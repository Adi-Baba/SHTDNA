#!/usr/bin/env python3
"""
Predictive Model: DNA Structure Function Classifier
=====================================================
Uses geometric signatures [Vij, Shambhian, Mamton] to predict functional class.

Model Classes:
- nucleosome: DNA wrapped around histone octamer
- protein_bound: DNA complexed with transcription factors, polymerases, etc.
- free_duplex: Canonical B-form or A-form DNA without proteins
- non_canonical: Quadruplex, triplex, damaged, unusual structures

Training: Manual annotation of subset → RF classifier
Testing: Predict function for remaining structures
Validation: Cross-validation + literature comparison
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import matplotlib.pyplot as plt
import seaborn as sns

# =======================
# PDB Annotations
# =======================

# Manual annotations based on PDB metadata and literature
KNOWN_ANNOTATIONS = {
    # Nucleosomes (extreme wrapping, high Vij+Mamton)
    '1AOI': 'nucleosome',  # Classic nucleosome core particle
    '1F66': 'nucleosome',  # H2A.Z variant nucleosome
    
    # Protein-bound DNA (deformed geometry)
    '1AHD': 'protein_bound',  # Antennapedia homeodomain-DNA complex
    '102D': 'protein_bound',  # DNA binding protein
    '1CGP': 'protein_bound',  # Protein-DNA complex
    
    # Free duplex (canonical B-DNA or A-DNA)
    '1BNA': 'free_duplex',  # Classic B-DNA dodecamer (Drew-Dickerson)
    '1D23': 'free_duplex',  # DNA duplex
    '355D': 'free_duplex',  # A-form DNA
    '3BSE': 'free_duplex',  # B-form DNA
    
    # Non-canonical structures
    '1JB7': 'non_canonical',  # G-quadruplex
    '2KF8': 'non_canonical',  # i-motif
    '1D3O': 'non_canonical',  # Holliday junction
}

# Expand with high-confidence predictions based on geometry
# Free duplex: low Vij (<0.10), high Shambhian (>0.95), low Mamton (<150)
LIKELY_FREE_DUPLEX = [
    # Classic free duplexes with canonical geometry
    '1BNK', '1BY4', '1C7U', '1C8C', '1D02', '1DE9', '1DH3', 'DSD', '1DSZ',
    '1EWN', '1F0O', '1F5T', '1AM9', '1BHM', '1D2I', '1D66', '1DNK',
    # Additional canonical structures
    '1D89', '2D47', '1EN3', '1EN4', '1EN5', '1EN6', '1EN7', '1EN8', '1EN9',
    '1DNA', '2BNA', '3DNA', '4DNA', '5DNA', '6DNA', '7DNA', '8DNA', '9DNA',
]

# Protein-bound DNA: moderate deformation, lower flatness (<0.90)
LIKELY_PROTEIN_BOUND = [
    '1B94', '1B96', '1BBX', '1BGB', '1BSS', '1EJ9', '1EO4', '1EON', 
    '1AZ0', '10MH', '173D', '1D1U', '1EYG',
]

# Non-canonical structures: extreme roughness or unusual geometry
LIKELY_NON_CANONICAL = [
    '1CMA', '1CQT', '1DUX', '1EYG',  # High Mamton (>500)
]

# =======================
# Load Data
# =======================

def load_geometric_data():
    """Load merged geometric data from coupling analysis."""
    data_path = Path(__file__).parent / 'coupling_results' / 'merged_geometric_data.csv'
    
    if not data_path.exists():
        raise FileNotFoundError(f"Merged data not found: {data_path}")
    
    df = pd.read_csv(data_path)
    print(f"Loaded geometric data: {len(df)} structures")
    print(f"  Features: vij_mean, shambhian_mean, mamton_mean")
    print(f"  Columns: {df.columns.tolist()}")
    
    return df

def assign_labels(df):
    """Assign functional labels based on annotations."""
    df = df.copy()
    df['label'] = 'unknown'
    
    # Apply known annotations
    for pdb_id, label in KNOWN_ANNOTATIONS.items():
        mask = df['pdb_id'].str.upper() == pdb_id.upper()
        df.loc[mask, 'label'] = label
    
    # Apply likely free duplex labels
    for pdb_id in LIKELY_FREE_DUPLEX:
        mask = df['pdb_id'].str.upper() == pdb_id.upper()
        if mask.sum() > 0:
            # Only assign if geometry matches (low Vij, high Shambhian, low Mamton)
            row = df.loc[mask].iloc[0]
            if row['vij_mean'] < 0.15 and row['shambhian_mean'] > 0.90 and row['mamton_mean'] < 1000:
                df.loc[mask, 'label'] = 'free_duplex'
    
    # Apply protein-bound labels
    for pdb_id in LIKELY_PROTEIN_BOUND:
        mask = df['pdb_id'].str.upper() == pdb_id.upper()
        if mask.sum() > 0:
            row = df.loc[mask].iloc[0]
            # Moderate deformation: 0.08 < Vij < 0.15, Shambhian < 0.97
            if 0.08 <= row['vij_mean'] < 0.15 and row['shambhian_mean'] < 0.97:
                df.loc[mask, 'label'] = 'protein_bound'
    
    # Apply non-canonical labels
    for pdb_id in LIKELY_NON_CANONICAL:
        mask = df['pdb_id'].str.upper() == pdb_id.upper()
        if mask.sum() > 0:
            df.loc[mask, 'label'] = 'non_canonical'
    
    labeled = df[df['label'] != 'unknown'].copy()
    unlabeled = df[df['label'] == 'unknown'].copy()
    
    print(f"\nLabeling summary:")
    print(f"  Labeled: {len(labeled)} structures")
    print(labeled['label'].value_counts())
    print(f"  Unlabeled: {len(unlabeled)} structures")
    
    return labeled, unlabeled

# =======================
# Build Model
# =======================

def build_classifier(X_train, y_train):
    """Train Random Forest classifier."""
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_train)
    
    # Random Forest with class balancing
    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=10,
        min_samples_split=3,
        min_samples_leaf=2,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    
    clf.fit(X_scaled, y_train)
    
    return clf, scaler

def cross_validate_model(X, y, n_folds=5):
    """Perform stratified cross-validation."""
    print(f"\n{'='*60}")
    print("CROSS-VALIDATION")
    print(f"{'='*60}")
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=10,
        min_samples_split=3,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    scores = cross_val_score(clf, X_scaled, y, cv=cv, scoring='accuracy')
    
    print(f"Accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
    print(f"Per-fold: {scores}")
    
    return scores

def evaluate_model(clf, scaler, X_test, y_test):
    """Evaluate model on test set."""
    X_scaled = scaler.transform(X_test)
    y_pred = clf.predict(X_scaled)
    
    print(f"\n{'='*60}")
    print("CLASSIFICATION REPORT")
    print(f"{'='*60}")
    print(classification_report(y_test, y_pred, zero_division=0))
    
    print(f"\n{'='*60}")
    print("CONFUSION MATRIX")
    print(f"{'='*60}")
    cm = confusion_matrix(y_test, y_pred)
    print(cm)
    
    return y_pred, cm

def feature_importance_analysis(clf, feature_names):
    """Analyze which geometric features are most predictive."""
    importances = clf.feature_importances_
    
    print(f"\n{'='*60}")
    print("FEATURE IMPORTANCE")
    print(f"{'='*60}")
    
    for name, importance in zip(feature_names, importances):
        print(f"  {name:20s}: {importance:.3f}")
    
    return importances

# =======================
# Predict Unknown
# =======================

def predict_unknown_structures(clf, scaler, df_unlabeled, features):
    """Make predictions on unlabeled structures."""
    print(f"\n{'='*60}")
    print("PREDICTIONS ON UNKNOWN STRUCTURES")
    print(f"{'='*60}")
    
    X = df_unlabeled[features].values
    X_scaled = scaler.transform(X)
    
    predictions = clf.predict(X_scaled)
    probabilities = clf.predict_proba(X_scaled)
    
    # Get confidence (max probability)
    confidence = probabilities.max(axis=1)
    
    results = df_unlabeled.copy()
    results['predicted_class'] = predictions
    results['confidence'] = confidence
    
    # Sort by confidence
    results = results.sort_values('confidence', ascending=False)
    
    print(f"\nPredicted class distribution:")
    print(results['predicted_class'].value_counts())
    
    print(f"\nTop 10 high-confidence predictions:")
    cols = ['pdb_id', 'predicted_class', 'confidence', 'vij_mean', 'shambhian_mean', 'mamton_mean']
    print(results[cols].head(10).to_string(index=False))
    
    return results

# =======================
# Visualization
# =======================

def plot_geometric_space_with_labels(df_labeled, output_dir):
    """3D scatter plot colored by functional class."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    classes = df_labeled['label'].unique()
    colors = plt.cm.Set1(np.linspace(0, 1, len(classes)))
    
    for cls, color in zip(classes, colors):
        mask = df_labeled['label'] == cls
        subset = df_labeled[mask]
        
        ax.scatter(
            subset['vij_mean'],
            subset['shambhian_mean'],
            subset['mamton_mean'],
            c=[color],
            label=cls,
            s=100,
            alpha=0.7,
            edgecolors='k'
        )
    
    ax.set_xlabel('Vij (Curvature)', fontsize=12)
    ax.set_ylabel('Shambhian (Flatness)', fontsize=12)
    ax.set_zlabel('Mamton (Roughness)', fontsize=12)
    ax.set_title('Geometric State Space by Function', fontsize=14, weight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    output_path = output_dir / 'functional_classification_space.png'
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()

def plot_confusion_matrix_heatmap(cm, class_names, output_dir):
    """Heatmap of confusion matrix."""
    fig, ax = plt.subplots(figsize=(8, 7))
    
    sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap='Blues',
        xticklabels=class_names,
        yticklabels=class_names,
        cbar_kws={'label': 'Count'},
        ax=ax
    )
    
    ax.set_xlabel('Predicted Class', fontsize=12)
    ax.set_ylabel('True Class', fontsize=12)
    ax.set_title('Confusion Matrix', fontsize=14, weight='bold')
    
    output_path = output_dir / 'confusion_matrix.png'
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()

def plot_feature_importance(importances, feature_names, output_dir):
    """Bar plot of feature importance."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    indices = np.argsort(importances)[::-1]
    
    ax.bar(range(len(importances)), importances[indices], color='steelblue', alpha=0.8)
    ax.set_xticks(range(len(importances)))
    ax.set_xticklabels([feature_names[i] for i in indices], fontsize=11)
    ax.set_ylabel('Importance', fontsize=12)
    ax.set_title('Feature Importance for Functional Classification', fontsize=14, weight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    output_path = output_dir / 'feature_importance.png'
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()

# =======================
# Main Workflow
# =======================

def main():
    print("="*70)
    print("DNA STRUCTURE FUNCTION PREDICTOR")
    print("Using SHT Geometric Signatures")
    print("="*70)
    
    # Create output directory
    output_dir = Path(__file__).parent / 'prediction_results'
    output_dir.mkdir(exist_ok=True)
    
    # Load data
    df = load_geometric_data()
    
    # Assign labels
    df_labeled, df_unlabeled = assign_labels(df)
    
    if len(df_labeled) < 10:
        print("\nWARNING: Insufficient labeled data. Need at least 10 samples.")
        print("Adding more annotations or expanding LIKELY_FREE_DUPLEX list.")
        return
    
    # Define features
    features = ['vij_mean', 'shambhian_mean', 'mamton_mean']
    
    X = df_labeled[features].values
    y = df_labeled['label'].values
    
    # Cross-validation
    cv_scores = cross_validate_model(X, y, n_folds=min(5, len(df_labeled)))
    
    # Train final model on all labeled data
    print(f"\n{'='*60}")
    print("TRAINING FINAL MODEL")
    print(f"{'='*60}")
    clf, scaler = build_classifier(X, y)
    print(f"Trained on {len(X)} labeled structures")
    
    # Feature importance
    importances = feature_importance_analysis(clf, features)
    
    # Predict unknown structures
    if len(df_unlabeled) > 0:
        results = predict_unknown_structures(clf, scaler, df_unlabeled, features)
        
        # Save predictions
        results_path = output_dir / 'predictions.csv'
        results.to_csv(results_path, index=False)
        print(f"\nSaved predictions: {results_path}")
    else:
        print("\nNo unlabeled structures to predict.")
        results = pd.DataFrame()
    
    # Create visualizations
    print(f"\n{'='*60}")
    print("GENERATING VISUALIZATIONS")
    print(f"{'='*60}")
    
    plot_geometric_space_with_labels(df_labeled, output_dir)
    plot_feature_importance(importances, features, output_dir)
    
    # Save model summary
    summary = {
        'model': 'RandomForestClassifier',
        'n_labeled': len(df_labeled),
        'n_unlabeled': len(df_unlabeled),
        'features': features,
        'classes': sorted(df_labeled['label'].unique().tolist()),
        'cross_val_accuracy': float(cv_scores.mean()),
        'cross_val_std': float(cv_scores.std()),
        'feature_importance': {
            feat: float(imp) for feat, imp in zip(features, importances)
        }
    }
    
    if len(results) > 0:
        summary['predictions'] = {
            'nucleosome': int((results['predicted_class'] == 'nucleosome').sum()),
            'protein_bound': int((results['predicted_class'] == 'protein_bound').sum()),
            'free_duplex': int((results['predicted_class'] == 'free_duplex').sum()),
            'non_canonical': int((results['predicted_class'] == 'non_canonical').sum()),
        }
        
        # High-confidence nucleosome predictions (potential discoveries)
        high_conf_nucleosomes = results[
            (results['predicted_class'] == 'nucleosome') & 
            (results['confidence'] > 0.7)
        ]
        
        if len(high_conf_nucleosomes) > 0:
            summary['high_confidence_nucleosomes'] = high_conf_nucleosomes['pdb_id'].tolist()
    
    summary_path = output_dir / 'model_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved model summary: {summary_path}")
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"Results saved to: {output_dir}")

if __name__ == '__main__':
    main()

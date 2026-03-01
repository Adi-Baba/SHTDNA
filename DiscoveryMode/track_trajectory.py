
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys

# Import SHT Analysis
sys.path.append(r"d:\OnlyHST\SHTDNA\PDB\code")
try:
    from sht_dna_analysis import compute_sigma_sht_from_pdb
except ImportError:
    # Fallback if relative import fails, try direct loading or assume function exists
    print("Could not import sht_dna_analysis. SHT extraction might fail.")

def analyze_trajectories():
    cloud_file = r"d:\OnlyHST\SHTDNA\DiscoveryMode\mechanical_state_cloud.csv"
    output_dir = Path(r"d:\OnlyHST\SHTDNA\DiscoveryMode")
    output_img = output_dir / "trajectory_map.png"
    
    # 1. Build the Manifold (PCA) from the Cloud
    print("Building Geometric Manifold...")
    df = pd.read_csv(cloud_file)
    features = ['v0_Curvature', 'Eta_r_Shambhian', 'Omega_Mamton']
    X_cloud = df[features].dropna().values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_cloud)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)
    
    # 2. Extract Walkers
    walkers = []
    
    # Walker 1: Ideal B-DNA
    path_ideal = r"d:\OnlyHST\SHTDNA\PDB\code\generated_pdbs\ideal_b_dna.pdb"
    walkers.append({'name': 'Ideal B-DNA', 'path': path_ideal, 'color': 'red', 'marker': 'o'})
    
    # Walker 2: 1BNA (Crystal)
    path_cryst = r"d:\OnlyHST\SHTDNA\PDB\dataset\hard_dataset\structures\1BNA.cif" 
    # Use standard 1BNA if available, otherwise any simple crystal
    walkers.append({'name': 'Crystal (1BNA)', 'path': path_cryst, 'color': 'cyan', 'marker': '^'})
    
    walker_results = []
    
    for w in walkers:
        print(f"Tracking Walker: {w['name']}...")
        try:
            # Run SHT Extraction
            res = compute_sigma_sht_from_pdb(w['path'], num_slices=20) # 20 slices for small DNA
            
            # Construct State Vector
            # v0 = |d^2r/ds^2| -> approximation from 'C_list' or usually we can't get v0 easily without the full pipeline.
            # Wait, extract_full_signature had logic for v0.
            # compute_sigma_sht_from_pdb returns (res_dict)
            # res_dict keys: 'C_axis', 'Q', 'S', etc.
            
            C = res['C_axis']
            S = res['S']
            Q = res['Q']
            
            # Calculate v0 (Curvature)
            # v0[i] = | C[i+1] - 2C[i] + C[i-1] | / (h^2) roughly
            # Let's do finite difference
            dC = np.gradient(C, axis=0) # Velocity
            ddC = np.gradient(dC, axis=0) # Acceleration (Curvature vector)
            v0 = np.linalg.norm(ddC, axis=1)
            
            # Calculate Eta_r (Flatness)
            # Eta = 1 - epsilon. Epsilon = sqrt( (l1-l2)^2 ... ) / (l1+l2)
            # Actually sht_dna_analysis computes 'eigenvalues' of Q.
            # lambda1, lambda2. ecc = (l1-l2)/(l1+l2). Eta = 1 - ecc.
            # Let's re-compute from Q
            eta_list = []
            omega_list = []
            
            for i in range(len(Q)):
                # Q[i] is [2 strands, 2, 2]
                # We average strands? Or sum?
                # Usually we sum the moments of the pair for the 'base pair' geometry
                # Q_bp = Q[i,0] + Q[i,1]
                Q_bp = np.sum(Q[i], axis=0)
                evals = np.linalg.eigvalsh(Q_bp)
                l1, l2 = evals[0], evals[1] # Sorted? eigvalsh sorts ascending
                # l1 < l2 usually?
                # eccentricity e = (l2 - l1) / (l2 + l1)
                # Eta = 1 - e = 1 - (l2-l1)/(l2+l1) = (2 l1) / (l1+l2) ??
                # Let's stick to definition: Eta = 1 - ecc
                if (l1+l2) > 1e-9:
                    ecc = abs(l1 - l2) / (l1 + l2)
                    eta = 1.0 - ecc
                else:
                    eta = 0.0
                eta_list.append(eta)
                
                # Omega (Mamton) from Skewness S
                # S[i] is [2 strands, 2] -> (Sx, Sy)
                # Sum strands: S_bp = S[i,0] + S[i,1]
                S_bp = np.sum(S[i], axis=0)
                omega = np.linalg.norm(S_bp)
                omega_list.append(omega)
            
            eta_arr = np.array(eta_list)
            omega_arr = np.array(omega_list)
            
            # Align lengths (gradient trims? no, purely same shape)
            
            # Form Data Frame
            w_df = pd.DataFrame({
                'v0_Curvature': v0,
                'Eta_r_Shambhian': eta_arr,
                'Omega_Mamton': omega_arr
            })
            
            # Project
            w_X = w_df[features].values
            w_X_scaled = scaler.transform(w_X)
            w_pca = pca.transform(w_X_scaled)
            
            walker_results.append({
                'name': w['name'],
                'data': w_pca,
                'color': w['color'],
                'marker': w['marker']
            })
            
        except Exception as e:
            print(f"Failed to track {w['name']}: {e}")

    # 3. Visualization
    plt.figure(figsize=(10, 8))
    
    # Background Density
    plt.hexbin(X_pca[:, 0], X_pca[:, 1], gridsize=50, cmap='gray_r', mincnt=1, alpha=0.3, label='Cloud Density')
    
    # Plot Walkers
    for w in walker_results:
        d = w['data']
        plt.plot(d[:, 0], d[:, 1], color=w['color'], linewidth=2, label=w['name'] + ' Path')
        plt.scatter(d[:, 0], d[:, 1], color=w['color'], marker=w['marker'], s=30)
        # Mark Start/End
        plt.text(d[0, 0], d[0, 1], 'Start', color=w['color'], fontsize=8)
        plt.text(d[-1, 0], d[-1, 1], 'End', color=w['color'], fontsize=8)

    plt.xlabel('Geometric Mode 1 (Integrity)', fontsize=12)
    plt.ylabel('Geometric Mode 2 (Shape)', fontsize=12)
    plt.title('Tracking DNA Trajectories in Reduced State Space', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_img, dpi=150)
    print(f"Trajectory Map saved to {output_img}")

if __name__ == "__main__":
    analyze_trajectories()

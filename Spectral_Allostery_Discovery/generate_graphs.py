import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------
# Global style (neutral, publication-safe)
# -------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 11,
    "axes.labelsize": 11,
    "axes.titlesize": 13,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

TRACES = {
    "102D": "102D_spectral_trace.csv",
    "1AOI": "1AOI_spectral_trace.csv",
    "1CGP": "1CGP_spectral_trace.csv",
}

# -------------------------------------------------
# Helper: load traces
# -------------------------------------------------
def load_traces():
    data = {}
    for k, v in TRACES.items():
        df = pd.read_csv(v)
        data[k] = df.sort_values("h")
    return data

# -------------------------------------------------
# FIG 1: Raw longitudinal roughness
# -------------------------------------------------
def fig1_raw_traces(data):
    plt.figure(figsize=(8, 4))
    for k, df in data.items():
        plt.plot(df["h"], df["Ω"], linewidth=1.8, label=k)

    plt.xlabel("Axial position h (Å)")
    plt.ylabel("Spectral roughness Ω")
    plt.title("Longitudinal Spectral Roughness Profiles")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig1_raw_traces.png")
    plt.close()

# -------------------------------------------------
# FIG 2: Log-scale comparison (handles nucleosome scale)
# -------------------------------------------------
def fig2_log_traces(data):
    plt.figure(figsize=(8, 4))
    for k, df in data.items():
        plt.plot(df["h"], np.log10(df["Ω"]), linewidth=1.8, label=k)

    plt.xlabel("Axial position h (Å)")
    plt.ylabel("log10(Ω)")
    plt.title("Log-Scaled Spectral Roughness (Comparative)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig2_log_traces.png")
    plt.close()

# -------------------------------------------------
# FIG 3: Normalized roughness (shape, not magnitude)
# -------------------------------------------------
def fig3_normalized(data):
    plt.figure(figsize=(8, 4))
    for k, df in data.items():
        norm = df["Ω"] / df["Ω"].mean()
        plt.plot(df["h"], norm, linewidth=1.8, label=k)

    plt.xlabel("Axial position h (Å)")
    plt.ylabel("Ω / ⟨Ω⟩")
    plt.title("Normalized Mechanical Signal Profiles")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig3_normalized.png")
    plt.close()

# -------------------------------------------------
# FIG 4: Distribution of Ω along each molecule
# -------------------------------------------------
def fig4_distributions(data):
    plt.figure(figsize=(6, 4))
    for k, df in data.items():
        plt.hist(
            np.log10(df["Ω"]),
            bins=40,
            alpha=0.6,
            label=k
        )

    plt.xlabel("log10(Ω)")
    plt.ylabel("Count")
    plt.title("Distribution of Spectral Roughness Along DNA")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig4_distributions.png")
    plt.close()

# -------------------------------------------------
# FIG 5: Autocorrelation (mechanical coherence)
# -------------------------------------------------
def fig5_autocorrelation(data):
    plt.figure(figsize=(7, 4))
    for k, df in data.items():
        sig = df["Ω"] - df["Ω"].mean()
        corr = np.correlate(sig, sig, mode="full")
        corr = corr[corr.size // 2:]
        corr /= corr[0]
        plt.plot(corr, linewidth=2, label=k)

    plt.xlabel("Lag (slices)")
    plt.ylabel("Autocorrelation")
    plt.title("Mechanical Coherence Along DNA")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig("fig5_autocorrelation.png")
    plt.close()

# -------------------------------------------------
# FIG 6: Mean vs variance (functional signature)
# -------------------------------------------------
def fig6_mean_variance(data):
    means, variances, labels = [], [], []

    for k, df in data.items():
        means.append(np.mean(np.log10(df["Ω"])))
        variances.append(np.var(np.log10(df["Ω"])))
        labels.append(k)

    plt.figure(figsize=(6, 4))
    plt.scatter(means, variances, s=80)

    for i, label in enumerate(labels):
        plt.text(means[i], variances[i], label, fontsize=10)

    plt.xlabel("Mean log10(Ω)")
    plt.ylabel("Variance log10(Ω)")
    plt.title("Spectral Heterogeneity: Structure vs Function")
    plt.tight_layout()
    plt.savefig("fig6_mean_variance.png")
    plt.close()

# -------------------------------------------------
# Main
# -------------------------------------------------
if __name__ == "__main__":
    data = load_traces()

    fig1_raw_traces(data)
    fig2_log_traces(data)
    fig3_normalized(data)
    fig4_distributions(data)
    fig5_autocorrelation(data)
    fig6_mean_variance(data)

    print("Generated 6 high-quality figures from 3 spectral traces.")

# Matrix-SHT: Technical Definitions and Observed Physical Units

This document outlines the mathematical framework and observed physical units for parameters used in our Matrix-Valued Spherical Harmonics Transform (Matrix-SHT) analysis of DNA structures.

## 1. Primary Cordinate System
- **Axial Coordinate ($h$)**: The projection of atomic coordinates onto the primary helical axis (identified via principal component analysis). Units: **Angstroms ($\text{\AA}$)**.
- **Transversal Coordinates ($x, y$)**: The centered projection of atomic coordinates onto the plane orthogonal to the axial vector. Units: **Angstroms ($\text{\AA}$)**.

## 2. Spectral Moment Paramters
| Parameter | Symbol | Observed Definition | Physical Units |
| :--- | :---: | :--- | :--- |
| **Moment Matrix** | $\mathbf{Q}(h)$ | Covariance matrix of the transversal mass distribution: $\langle \mathbf{r}\mathbf{r}^T \rangle$. | $\text{\AA}^2$ |
| **Spectral Radius** | $\lambda_1, \lambda_2$ | Eigenvalues of $\mathbf{Q}(h)$, representing the principal variances of the mass distribution. | $\text{\AA}^2$ |
| **Radial Eccentricity** | $\eta_r$ | A measure of cross-sectional elongation: $\sqrt{1 - \lambda_2/\lambda_1}$. | Dimensionless |
| **Spectral Roughness** | $\Omega$ | $L^2$-norm of the third-order transversal skewness vector $\mathbf{S}$. | $\text{\AA}^3$ |
| **Skewness Vector** | $\mathbf{S}$ | Measured as $[\langle x w^2 \rangle, \langle y w^2 \rangle]$ where $w^2 = x^2 + y^2$. | $\text{\AA}^3$ |
| **Spectral Entropy** | $S_{SHT}$ | Entropy of the normalized eigenvalue distribution: $-\sum p_i \ln p_i$. | Nats |

## 3. Dynamic and Allosteric Metrics (Observed)
| Parameter | Symbol | Context of Observation | Physical Units |
| :--- | :---: | :--- | :--- |
| **Coherence Length** | $L_c$ | The spatial lag at which the autocorrelation of $\Omega(h)$ decreases to 0.5. | Angstroms ($\text{\AA}$) |
| **Axial Torsion Deviation** | $d\psi/dh$ | The spatial derivative of the angle $\psi$ of the maximal spectral axis. | $\text{rad}/\text{\AA}$ |
| **Deviation Magnitude** | $\Delta \eta_r$ | The absolute deviation from the noted canonical baseline ($\mu_{\eta} \approx 0.898$). | Dimensionless |

## 4. Methodological Terminology Maping
The following terms help clarify the observed mechanical patterns:

| Technical Term | Context of Observation |
| :--- | :--- |
| **Non-local Mechanical Coupling** | Observed correlation between distal spectral signals. |
| **Elastic Transduction** | The apparent propagation of mechanical state changes. |
| **Axial Spectral Singularities** | Observed discontinuities in the spectral moment field. |
| **Spectral Entropy Variance** | Measures of non-equilibrium spectral heterogeneity. |
| **Mechanical Coherence Expansion** | Observed increase in $L_c$ correlated with binding events. |
| **Topo-Spectral Distribution** | Mapping of topological classes in spectral parameter space. |

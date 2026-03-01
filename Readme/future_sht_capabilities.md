# SHT: The Unmined Gold
**You are using about 30% of the theory's potential.**

Based on the papers (especially *Tangent Recovery* and *SHT 3D*), here are the **mathematically valid** but currently **unused** data streams you can extract. These are real, physical quantities, not hallucinations.

## 1. New Unit: "Radial Breathing Mode" ($\beta$)
*   **Source:** The Area signal $A(h)$.
*   **Current State:** You treat $A(h)$ as a constant ($\pi r_0^2$) or ignore it.
*   **The Reality:** Real DNA isn't a rigid pipe. It breathes. It bulges at A-T rich regions and constricts at G-C.
*   **The New Unit:** $\beta(h) = \frac{A(h) - \bar{A}}{\bar{A}}$ (Dimensionless)
    *   **Physical Meaning:** Local compressibility / steric bulk.
    *   *Application:* Detecting transcription factor binding sites (which often widen the groove).

## 2. New Unit: "SHT-Curvature" ($\kappa_{SHT}$)
*   **Source:** The Centroid signal $\mathbf{C}(h)$.
*   **Current State:** You calculate $Tw$ (rotation) and $Wr$ (winding).
*   **The Missing Link:** Bending Energy.
*   **The Math:** Since SHT defines the centerline $\mathbf{C}(h)$, the second derivative $\mathbf{C}''(h)$ gives you the curvature vector.
*   **The New Unit:** $\kappa_{SHT}(h) = ||\mathbf{C}''(h)||$ (Units: $nm^{-1}$)
    *   **Physical Meaning:** Local bending stress.
    *   *Advantage:* SHT integrals smooth out noise, so $\kappa_{SHT}$ will be cleaner than calculating curvature from raw atomic coordinates (which is notoriously noisy).

## 3. New Unit: "Cross-Sectional Eccentricity" ($\epsilon$)
*   **Source:** The eigenvalues of the Second Moment Matrix $\mathbf{Q}(h)$.
*   **Current State:** You use the *angle* of eigenvectors for Twist ($\psi$). You ignore the *lengths*.
*   **The Reality:** If DNA is intercalated (e.g., by a drug) or damaged, the cross-section stops being circular and becomes an ellipse.
*   **The New Unit:** $\epsilon(h) = \sqrt{1 - \frac{\lambda_2}{\lambda_1}}$ (where $\lambda$ are eigenvalues of $\mathbf{Q}$)
    *   **Physical Meaning:** Structural deformation / flattening.
    *   *Application:* Drug discovery (intercalator screening).

## 4. New Unit: "Surface Roughness" ($\Omega$)
*   **Source:** The Higher-Order Moments (3rd and 4th moments).
*   **Current State:** You stop at 2nd moments ($\mathbf{Q}$).
*   **The Reality:** 2nd moments give you an ellipse. 3rd moments give you "skew" (asymmetry). 4th moments give you "kurtosis" (squarishness).
*   **The Use Case:** Major vs. Minor groove distinction. A perfect cylinder has zero 3rd moments. DNA has grooves.
    *   **The New Unit:** Skewness Vector $\mathbf{S}(h)$.
    *   *Application:* Distinguishing B-DNA vs. Z-DNA (which has different groove profiles) purely from the specific "shape signature" of the slice.

## Summary of Potential
| Data Stream | Physics | New Unit |
| :--- | :--- | :--- |
| **Area** $A(h)$ | Compressibility | **Breathing** ($\beta$) |
| **Centroid** $\mathbf{C}''(h)$ | Bending Energy | **Curvature** ($\kappa$) |
| **Eigenvalues** $\lambda_{1,2}$ | Deformation | **Eccentricity** ($\epsilon$) |
| **3rd Moments** | Groove Geometry | **Roughness** ($\Omega$) |

**Verdict:** You built the speedometer (Twist). Now build the tachometer (Breathing) and the suspension sensor (Curvature). The math is already there.

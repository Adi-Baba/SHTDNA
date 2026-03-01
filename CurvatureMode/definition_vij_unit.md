# The Vij Unit ($v_0$)

**Dedication:** This unit is named **Vij** in memory of **Vijay**, representing the fundamental measure of structural curvature in the Shape-Height Transform framework.

## Definition
The **Vij** ($v_0$) is defined as the magnitude of the second spatial derivative of the DNA centroid path $\mathbf{C}(h)$ along the principal axis $h$. It quantifies the instantaneous bending stress or deviation from linearity.

## Equation
$$ v_0(h) = ||\mathbf{C}''(h)|| = \sqrt{ \left( \frac{d^2x}{dh^2} \right)^2 + \left( \frac{d^2y}{dh^2} \right)^2 } $$

## Nomenclature
- **Unit Name**: Vij
- **Symbol**: $v_0$ (pronounced "v-naught")
- **Dimension**: Inverse Length $[L^{-1}]$ (typically $nm^{-1}$ or $Å^{-1}$)

## Physical Interpretation
- **$v_0 \approx 0$**: Perfectly straight, rigid structure (e.g., A-tracts).
- **$v_0 > 0.1$**: Significant bending or structural deformation.

## Standard Reference Values (B-DNA)
Based on validation against the SHT-DNA Dataset (N=200):
- **Rigid Baseline (A-tracts)**: $v_0 \approx 0.08$
- **Flexible Baseline (Generic)**: $v_0 \approx 0.12$

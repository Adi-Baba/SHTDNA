
import numpy as np
from scipy.optimize import curve_fit

# Data from lost_stiffness_log.txt
# Midpoints of the threshold bins? No, the log uses "greater than" thresholds. 
# Let's treat the threshold as the characteristic 'Eta_r' of that regime.
eta = np.array([0.700, 0.715, 0.729, 0.744, 0.759, 0.774, 0.788, 
                0.803, 0.818, 0.833, 0.847, 0.862, 0.877, 0.892, 
                0.906, 0.921, 0.936, 0.951, 0.965, 0.980])

Lp = np.array([7.23, 7.56, 8.11, 8.61, 9.25, 9.94, 10.98,
               12.43, 13.70, 14.78, 15.96, 16.78, 18.07, 20.29,
               22.75, 26.46, 37.69, 58.48, 393.05, 657.44])

# Fit log(Lp) to handle the 7 -> 600 range

# Model 1: Exponential
# ln(Lp) = ln(A) + k * x
def model_exp_log(x, ln_a, k):
    return ln_a + k * x

# Model 2: Critical Divergence
# ln(Lp) = ln(A) - gamma * ln(1 - x)
def model_crit_log(x, ln_a, gamma):
    return ln_a - gamma * np.log(1 - x + 1e-9)

print("Fitting Log-Scaling Laws...")
try:
    log_Lp = np.log(Lp)
    
    # Fit Exp
    popt_exp, _ = curve_fit(model_exp_log, eta, log_Lp, p0=[0, 5])
    ss_exp = np.sum((model_exp_log(eta, *popt_exp) - log_Lp)**2)
    A_exp = np.exp(popt_exp[0])
    
    print(f"Exponential: Lp = {A_exp:.4f} * exp({popt_exp[1]:.4f} * eta)")
    print(f"Log-SSR: {ss_exp:.4f}")
    
    # Fit Critical
    popt_crit, _ = curve_fit(model_crit_log, eta, log_Lp, p0=[2, 1])
    ss_crit = np.sum((model_crit_log(eta, *popt_crit) - log_Lp)**2)
    A_crit = np.exp(popt_crit[0])
    
    print(f"Critical:    Lp = {A_crit:.4f} * (1 - eta)^(-{popt_crit[1]:.4f})")
    print(f"Log-SSR: {ss_crit:.4f}")

    if ss_crit < ss_exp:
        print("\nWINNER: Critical Divergence (Phase Transition).")
    else:
        print("\nWINNER: Exponential Scaling.")
        
except Exception as e:
    print(e)

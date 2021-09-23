import numpy as np

def beta_cooperation(beta_max, coop_frac, type_response="exponential"):

    if type_response=="exponential":
        beta = beta_max*np.exp(-coop_frac)
    elif type_response=="linear":
        beta = beta_max*(1-coop_frac)  + np.exp(-1)*beta_max*coop_frac
    elif type_response=="concave":
        p = 1/(4*(beta_max*np.exp(-1) - beta_max))
        beta2 = (coop_frac**2)/(4*p) + beta_max
        beta  = beta2
    elif type_response=='convex':
        p = 1/(4*(beta_max - beta_max*np.exp(-1)))
        beta2 = ((coop_frac-1)**2)/(4*p) + beta_max*np.exp(-1)
        beta  = beta2
    elif type_response=='s-shape':
        beta_min = beta_max*np.exp(-1)
        r = 10
        c = beta_max / (1 + (beta_max/beta_min - 1)/np.exp(-r))
        beta3 = beta_max/(1 + ((beta_max - beta_min)/beta_min) * np.exp(-(1-coop_frac)*r))
        beta = beta3
    # elif type_response=='s-shape-inv':
        # beta_y = beta_max*(-coop_frac)  + np.exp(-1)*beta_max*coop_frac
        # beta_min = beta_max*np.exp(-1)
        # # r = 10
        # # c = beta_max / (1 + (beta_max/beta_min - 1)/np.exp(-r))
        # # beta3 = (np.log(beta_max/beta_y - 1) - np.log((beta_max - beta_min)/beta_min))/r +1
        # r = 10
        # c = beta_max / (1 + (beta_max/beta_min - 1)/np.exp(-r))
        # beta3 = beta_max/(1 + ((beta_max - beta_min)/beta_min) * np.exp((1-coop_frac)*r))
        # beta = beta3

    return beta
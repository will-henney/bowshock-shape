import numpy as np

def alpha_from_c_beta(c, beta):
    sb = np.sqrt(beta)
    return c/(1.0 + sb) + (1.0 + 2.0*sb)/6.0


def c_from_beta_k(beta, k):
    return (1.0 - beta - 3*3.0*k/4.0)/30.0


def k_from_xi(xi):
    return 2.0/xi - 2.0


def alpha_from_beta_xi(beta, xi):
    return alpha_from_c_beta(c_from_beta_k(beta, k_from_xi(xi)), beta)

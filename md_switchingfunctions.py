"""
A collection of switching functions used by different MD engines
"""

import numpy as np

def switch_CHARMM(r, r_in, r_out):
    """

    Function:
    ---------
    ``S(r) = (r_out^2 - r^2)^2 * (r_out^2 + 2*r^2 - 3*r_in^2) / (r_out^2 - r_in^2)^3``

    Parameters:
    -----------
    r : numpy.ndarray
        distance
    r_in : float
        inner cutoff
    r_out : float
        outer cutoff

    Returns:
    --------
    s : numpy.ndarray
        switching function
    """

    s = np.ones(r.size)

    mask = (r_in < r) & (r < r_out)
    s[mask] = np.power(r_out ** 2 - np.power(r[mask], 2), 2)
    s[mask] *= (r_out ** 2 + 2 * np.power(r[mask], 2) - 3 * r_in ** 2)
    s[mask] /= np.power(r_out ** 2 - r_in ** 2, 3)

    s[r >= r_out] = 0

    return s


def switch_openmm(r, r_switch, r_cutoff):
    """
    Function:
    ---------
    ``S(x) = 1 -6*x^5 + 15*x^4 - 10*x^3``
    ``x(r) = (r - r_switch)/(r_cutoff - r_switch)``

    Parameters:
    -----------
    r : numpy.ndarray
        distance
    r_switch : float
        inner cutoff
    r_cutoff : float
        outer cutoff
    Returns:
    --------
    s : numpy.ndarray
        switching function
    """
    s = np.ones(r.size)
    mask = (r_switch < r) & (r < r_cutoff)

    x = (r[mask] - r_switch) / (r_cutoff - r_switch)
    s[mask] *= 1 - 6 * np.power(x, 5) + 15 * np.power(x, 4) - 10 * np.power(x, 3)
    s[r >= r_cutoff] = 0
    return s


def switch_mdf(r, r_m, r_cut):
    """
    Function:
    ---------
    ``S(r) = (1 - x^3) * (1 + 3*x + 6*x^2)``
    ``x = (r - r_m) / (r_cut - r_m)

    Parameters:
    -----------
    r : numpy.ndarray
        distance
    r_m : float
        inner cutoff
    r_cut : float
        outer cutoff
    Returns:
    --------
    s : numpy.ndarray
        switching function
    """
    s = np.ones(r.size)
    mask = (r_m < r) & (r < r_cut)

    x = (r[mask] - r_m) / (r_cut - r_m)
    s[mask] *= np.power(1 - x, 3) * (1 + 3 * x + 6 * np.power(x, 2))
    s[r >= r_cut] = 0
    return s
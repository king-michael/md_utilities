#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
PLUMED switchingfunctions

see: http://plumed.github.io/doc-master/user-doc/html/switchingfunction.html

Author : Michael King <michael.king@uni-konstanz.de>
"""

import numpy as np

__author__ = "Michael King"


def rational(r,  # Type: Union[List, Tuple, numpy.ndarray]
             r_0,  # Type: float
             d_0=0.0,  # Type: float
             n=6,  # Type: int
             m=None,  # Type: Union[Int, None]
             d_max=None,  # Type: Union[Float, None]
             ):  # Type: numpy.ndarray
    """
    RATIONAL switching function from PLUMED

    `s(r) = ( 1 - ( (r-d_0) / r_0 )^n ) / ( 1 - ( (r-d_0) / r_0 )^m )`

    Parameters
    ----------
    r : list or tuple or numpy.ndarray
    r_0 : float
    d_0 : float, optional
    n : int, optional
    m : int or None, optional
        `None` == `m = 2n`
        (Default is `None`.)
    d_max : float or None
        (Default is `None`.)

    Return
    ------
    s : numpy.ndarray
        rational switchingfunction applied on `r`
    """
    r = np.asarray(r)

    if m is None:
        m = 2*n

    rdist = (r - d_0) / d_0

    s = (1 - np.power(rdist, n)) / (1 - np.power(rdist, m))

    if d_max is not None:
        s[r >= d_max] = 0.0

    return s

"""
Implementation of basic math function is python.
"""

from functools import reduce


def gcd(a, b):
    """
    Return greatest common divisor using Euclid's Algorithm.

    Parameters
    ----------
    a : int
        integer
    b : int
        integer

    Returns
    -------
    int
        greatest common divisor
    """
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """
    Return lowest common multiple.

    Parameters
    ----------
    a : int
        integer
    b : int
        integer

    Returns
    -------
    int
        lowest common multiple
    """
    return a * b // gcd(a, b)


def lcmm(*args):
    """
    Return lowest common multiple of Integers.

    Parameters
    ----------
    *args : Iterable[int]
        Integers

    Returns
    -------
    int
        lowest common multiple
    """
    """Return lcm of args."""
    return reduce(lcm, args)

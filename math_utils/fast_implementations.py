import numpy as np


def np_unique_int(array, return_counts=False):
    """
    Fast variant of ``np.unique(array, return_counts=True)``

    Only works with integer values.

    Parameter
    ---------
    array : np.ndarray
        Input array. Has to be 1-D.

    return_counts : bool, optional
        Return the counts next to the unique values.
        Default is `False`.

    Returns
    -------
    unique : np.ndarray
        Unique elements
    unique_counts : np.ndarray, optional
        Counts of Unique elements

    Raises
    ------
    TypeError
        If the input `array.dtype == float` raises an TypeError.
        This can be avoided by `array.astype(int)`.

    Examples
    --------
    >>> data = np.array([1,2,3,2,3])
    >>> np_unique_int(data)
    array([1, 2, 3])


    >>> data = np.array([1,2,3,2,3])
    >>> np_unique_int(data, return_counts=True)
    (array([1, 2, 3]), array([1, 2, 2]))

    >>> data = np.array([1,2,3,2,3], dtype=np.float64)
    >>> np_unique_int(data.astype(int))
    array([1, 2, 3])
    """

    bincount = np.bincount(array.astype(int))
    unique = np.where(bincount.astype(np.bool))[0]

    if return_counts:
        unique_counts = bincount[unique]
        return unique, unique_counts
    return unique

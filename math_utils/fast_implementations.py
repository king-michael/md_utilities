import numpy as np

def get_counts_ids(array):
    """
    Fast variant of ``np.unique(array, return_counts=True)``

    Only works with integer values.

    Parameter
    ---------
    array : np.ndarray

    Returns
    -------
    unique : np.ndarray
        Unique elements
    unique_counts : np.ndarray
        Counts of Unique elements
    """
    # ids, counts = np.unique(cluster_id2[-1], return_counts=True)
    bincount = np.bincount(array.astype(int))
    unique = np.where(bincount.astype(np.bool))[0]
    unique_counts = bincount[unique]
    return unique, unique_counts


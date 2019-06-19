import numpy as np
import MDAnalysis

def drop_duplicated_frames(u, data, shift=0, axis=0):
    """
    Drops duplicated frames from the data.
    Useful if multiple trajectories with overlapping frames are loaded.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Used MDAnalysis Universe to create the data
    data : np.ndarray
        Data
    shift : int
        Shift the index by the given amount. Default is ``0``.
    axis : int, optional
        Axis where to drop the data. Default is ``0``.

    Returns
    -------
    data : np.ndarray

    """
    times_duplicate = (np.cumsum(u.trajectory.total_times / u.trajectory.dts)-shift).astype(int)
    return np.delete(data, times_duplicate[times_duplicate < data.shape[axis]], axis=axis )


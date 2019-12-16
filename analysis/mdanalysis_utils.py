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


def center_around_group(ag, iterations=1, center='origin', weights='mass'):
    """
    Transformation to center around a group.

    Parameters
    ----------
    ag : MDAnalysis.core.groups.AtomGroup
        Atomgroup to center around
    iterations : int, optional
        Number of iterations used. Default is `1`.
    center : str, optional
        Where to place the center. `origin` or 'box'. Default is `origin`.
    weights : str or None or array, optional
        Weights to use. Default is `None`.
        If `weigths` is None center_of_geometry is used.
        If `weights` is `mass` center_of_mass is used.
        If `weights` is an array the array is used.


    Returns
    -------

    """
    from MDAnalysis.lib.distances import apply_PBC

    ix = ag._ix
    center = True if center is 'origin' else False

    if weights is None:
        weights = np.ones((ix.shape,1), dtype=np.float64)
    elif weights == 'mass':
        weights = ag.masses[:, None].astype(np.float64)
    else:
        weights = np.atleast_2d(weights.astype(np.float64))
        if weights.shape[0] == 1:
            weights = weights.T

    sum_weights = weights.sum()

    def transformation(ts):
        box_half = 0.5 * ts.dimensions[:3]

        for i in range(iterations):
            # wrap atoms back into box
            pos_center = apply_PBC(ts.positions[ix], ts.dimensions)
            # get center of mass
            com = np.sum(pos_center * weights, 0) / sum_weights
            # com = self.ag.center_of_mass(pbc=False)

            # center around the group
            ts.positions = apply_PBC(ts.positions, ts.dimensions)
            ts.positions -= com - box_half
            ts.positions = apply_PBC(ts.positions, ts.dimensions)
            # move center of atoms to origin
            #ts.positions -= box_half
        if center:
            ts.positions -= box_half
        return ts

    return transformation
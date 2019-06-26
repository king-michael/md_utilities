import numpy as np
from MDAnalysis.lib.distances import distance_array, self_distance_array


def calculate_cluster_per_atom(ag, cutoff, backend='serial'):
    """
    Distance based clustering.

    Faster numpy implementation, gives the same results as ``cluster/atom`` in LAMMPS .

    Speedup for 100 atoms:
        * this implementation : ``864 µs ± 7.63 µs``
        * original code (written as in C) : ``7.4 s ± 47.6 ms``

    Parameters
    ----------
    ag : MDAnalysis.core.groups.AtomGroup
        Atom group to cluster
    cutoff : float
        Cutoff distance
    backend : str, optional
        Which backend to use for the distance calculation.
        Note: ``serial`` may be faster then ``OpenMP`` for smaller systems.
        Default is ``serial``.

    Returns
    -------
    clusterID : np.ndarray
        Array with the cluster ID's of shape ``(ag.n_atoms,)``.
        The cluster ID is the lowest ID of the atoms in the cluster.
    """
    clusterID = ag.ids + 1

    n_atoms = ag.n_atoms
    distances = self_distance_array(ag.positions, box=ag.ts.dimensions, backend=backend)
    ids = np.where(distances < cutoff)[0]

    mi, mj = np.array(np.triu_indices(n_atoms,k=1))[:,ids]

    while True:
        done = True
        mask = np.where(clusterID[mi] != clusterID[mj])[0]
        if mask.size > 0:
            i, j = mi[mask], mj[mask]
            clusterID[i] = clusterID[j] = np.min(np.asarray([clusterID[i], clusterID[j]]), axis=0)
            done = False
        if done:
            break
    return clusterID


def _calculate_cluster_per_atom_lammps(ag, cutoff):
    """
    Distance based clustering.

    Same implementation as ``cluster/atom`` in LAMMPS.  (written as in C)

    For 100 atoms `7.4 s ± 47.6 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)`

    Parameters
    ----------
    ag : MDAnalysis.core.groups.AtomGroup
        Atom group to cluster
    cutoff : float
        Cutoff distance

    Returns
    -------
    clusterID : np.ndarray
        Array with the cluster ID's of shape ``(ag.n_atoms,)``.
        The cluster ID is the lowest ID of the atoms in the cluster.
    """
    clusterID = ag.ids + 1
    while True:
        done = True
        for i in range(0, ag.n_atoms):

            xtmp = ag[i].position

            for j in range(0, ag.n_atoms):
                if clusterID[i] == clusterID[j]:
                    continue

                r = distance_array(xtmp, ag[j].position, box=ag.ts.dimensions)[0][0]
                if r < cutoff:
                    clusterID[i] = clusterID[j] = min(clusterID[i], clusterID[j])
                    done = 0;
        if done:
            break
    return clusterID



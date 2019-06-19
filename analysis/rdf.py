import numpy as np
from MDAnalysis.lib.distances import self_distance_array
from progress_reporter import ProgressReporter_


def human_readable_to_bytes(size):

    """Given a human-readable byte string (e.g. 2G, 10GB, 30MB, 20KB),
      return the number of bytes.  Will return 0 if the argument has
      unexpected form.
   """
    if (size[-1] == 'B'):
        size = size[:-1]
    if (size.isdigit()):
        bytes = int(size)
    else:
        bytes = size[:-1]
        unit = size[-1]
        if (bytes.isdigit()):
            bytes = int(bytes)
            if (unit == 'G'):
                bytes *= 1073741824
            elif (unit == 'M'):
                bytes *= 1048576
            elif (unit == 'K'):
                bytes *= 1024
            else:
                bytes = 0
        else:
            bytes = 0
    return bytes ,size +'B'


def calculate_RDF(ag,
                  range=(0, 15),
                  bins=100,
                  start=0,
                  end=None,
                  step=1,
                  backend='serial',  # type: str {'serial', 'OpenMP'}
                  max_memory_usage=None,
                  verbose=False
                  ):
    """
    Calculates the radial distribution function for all atoms in the AtomGroup with them self.

    Parameter
    ---------
    ag : MDAnalysis.core.groups.AtomGroup
        Atom group
    range : tuple(int,int), optional
        Bin range to use. Default is `(0, 15)`.
    bins : int, optional
        Number of bins used. Default is `100`.
    start : int, optional
        Starting frame. Default is `0`.
    end : int or None, optional
        Final frame. `None` for last frame. Default is `None`.
    step : int, optional
        Step size. Default is `1`.
    backend : str, optional
        Backend to use.  `{'serial', 'OpenMP'}`. Default is `serial`.
    max_memory_usage : int or None, optional
        Maximum memory to use.
        If it's not `None`, results will be buffered before calculating the histogram.
        Default is `None`.

    Returns
    -------
    bins : numpy.ndarray
        Array of the bin centers.
    rdf : numpy.ndarray
        Radial distribution function.

    Examples
    --------
    >>> bins, rdf = calculate_RDF(ag, bins=100, range=(0, 15),
                                  start=0, end=1000,
                                  backend='OpenMP',
                                  max_memory_usage=2 * 1024 * 1024 * 1024, # 2 GB
                                  verbose=True)
    """

    # settings
    rdf_settings = dict(bins=bins,
                        range=range)

    # general constants
    n_frames = len(ag.universe.trajectory[start:end:step])
    n_pairs = ag.n_atoms * (ag.n_atoms - 1) // 2

    # handle maximum memory requirement
    if max_memory_usage is not None:
        # calculate the needed array size
        size_distance_array = int(n_pairs * 8)
        n_buffer = max_memory_usage // size_distance_array
        assert n_buffer != 0, "Not enough memory to calculate the rdf with this implementation. " + \
                              "Need at least {} bytes".format(size_distance_array)
    else:
        n_buffer = 1

    # init storage array
    dummy_storage = np.empty((n_buffer, n_pairs), dtype=np.float64)

    # init histogram
    count, edges = np.histogram([-1], **rdf_settings)

    # initialze the volume
    volume = 0

    # initialize buffer
    b = 0
    if verbose:
        p = ProgressReporter_()
        p.register(n_frames, description="calculate RDF")
    for ts in ag.universe.trajectory[start:end:step]:
        if verbose:
            p.update(1)
        # calculate distances
        self_distance_array(ag.positions, box=ts.dimensions,
                            result=dummy_storage[b], backend=backend)
        b += 1  # go to the next buffer
        if b == n_buffer:
            tmp_count, _ = np.histogram(dummy_storage, **rdf_settings)
            count += tmp_count
            b = 0  # reset b

        volume += ts.volume
    if b > 0:
        tmp_count, _ = np.histogram(dummy_storage[:b], **rdf_settings)
        count += tmp_count
    if verbose:
        p.finish()

    # Volume in each radial shell
    vol = np.power(edges[1:], 3) - np.power(edges[:-1], 3)
    vol *= 4 / 3.0 * np.pi

    # Average number density
    box_vol = volume / n_frames
    density = n_pairs / box_vol

    rdf = count / (density * vol * n_frames)
    bins = (edges[:-1] + edges[1:]) / 2.0

    return bins, rdf


def calculate_S(r, gr, q=None,  rho=1.0):
    r"""
    Calculate the isotropic structure factor

    `S(q) = 1 + (4*pi*rho)/q * sum( r * np.sin(r*q) * [g(r) - 1] )`

    Parameters
    ----------
    r : np.ndarray
        Radius
    gr : np.ndarray
        Radial distribution function.
    q : np.ndarray or None
        Vector of magintudes of wave vectors for which S(q) should be calculated.
        If `None` a wave vectors in the interval `[1.0/r.max(), 1.0/r.min())`
    rho : float
        density. For a isotropic system its a scalar.

    Returns
    -------
    q : np.ndarray
        Vector of magintudes of the wave vectors
    s : np.ndarray
        S(q)
    """
    # define wave vectors if not predefined.
    if q is None:
        q = np.linspace(1.0 / bins.max(), 1.0 / bins.min(), num=100)

    # define a grid
    qq, rr = np.meshgrid(q, r)
    _, grgr = np.meshgrid(q, gr)

    # calculate s
    s = 1.0 + ((4.0 * np.pi * rho) / q)
    # add the integral over r for every q
    s *= np.sum(rr * np.sin(qq * rr) * (grgr - 1), axis=0)

    return q, s
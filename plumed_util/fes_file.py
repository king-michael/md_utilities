"""
Function to load a Free Energy File produced by
``plumed sum_hills``

author : Michael King
"""

import numpy as np
from collections import namedtuple


def get_header(inp_file):
    # Type: (str) -> List[namedtuple]
    """
    Returns the header of a `fes.dat` file.

    Parameters
    ----------
    inp_file : str
        path to the ``fes.dat`` file

    Returns
    -------
    header : list[namedtuple]
        List of the CV's with:
        ``('name', 'min', 'max', 'nbin', 'periodic')``
    """
    header = []
    CV = namedtuple('CV', ['name', 'min', 'max', 'nbin', 'periodic'])
    with open(inp_file, 'r') as fp:
        firstline = fp.readline().split()
        if firstline[-1] == 'projection':
            NUM_CV = int((len(firstline) - 1) / 2)
        else:
            NUM_CV =int((len(firstline) - 3) / 2)
        handle_pi = lambda x: float(eval(x.replace('pi', '*pi') if re.match('\dpi', x) else x))
        for i in range(NUM_CV):
            header.append(CV(name=firstline[i+2],
                min=handle_pi(fp.readline().split()[-1]),
                max=handle_pi(fp.readline().split()[-1]),
                nbin=int((fp.readline()).split()[-1]),
                periodic=True if str(fp.readline().split()[-1]) == 'true' else False))
    return header


def load(inp_file):
    # Type: (str) -> Tuple[list, numpy.ndarray]
    """
    Function to load the ``fes.dat`` file produced by ``PLUMED``.


    Parameters
    ----------
    inp_file : str
        path to the ``fes.dat`` file

    Returns
    -------
    header : list[namedtuple]
        List of the CV's with:
        ``('name', 'min', 'max', 'nbin', 'periodic')``
    edges : list[numpy.ndarray]
        List of edges
    data : numpy.ndarray
        fes data
    """

    header = get_header(inp_file)
    num_cv = len(header)

    raw_data = np.genfromtxt(inp_file, usecols=[i for i in range(num_cv + 1)])
    edges = [raw_data[:, i].reshape([cv.nbin for cv in header[::-1]])[
                 tuple(slice(None) if j == i else 0 for j in range(num_cv - 1, -1, -1))]
             for i in range(num_cv)]

    data = raw_data[:, num_cv].reshape([cv.nbin for cv in header[::-1]])

    return header, edges, data.T


def plot(header, data, ax=None):
    # Type: (list, numpy.ndarray) -> matplotlib.image.AxesImage
    """

    Function to plot a ``fes.dat`` file.

    Parameters
    ----------
    header : list
        header of the ``fes.dat`` file
    data : numpy.ndarray
        fes.dat
    ax : None or matplotlib.axes._subplots.AxesSubplot
        ``ax`` object to plot in.
        If ``None`` an ax object will be created.
            * Get ``Figure`` with: ``fig = img.get_figure()``
            * Get ``Axes`` with: ``ax = fig.get_axes()``

    Returns
    -------
    img : matplotlib.image.AxesImage

    """

    import matplotlib.pyplot as plt
    from itertools import chain

    if ax is None:
        fig, ax = plt.subplots()

    extent = list(chain.from_iterable(map(lambda x: [x.min, x.max], header)))

    img = ax.imshow(data, origin='lower',
               extent=extent,
               aspect=float(extent[1] - extent[0]) / (extent[3] - extent[2])
               )

    return img
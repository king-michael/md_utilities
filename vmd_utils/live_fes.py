"""
Functions to plot a free energy surface in VMD and connect it with vmd_frame changes.

TODO:
 - interpolation for 1D surfaces. Currently it only takes the very next point.
"""

from VMD import evaltcl
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from collections import namedtuple
import os
import re


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
            NUM_CV = int((len(firstline) - 3) / 2)
        handle_pi = lambda x: float(eval(x.replace('pi', '*pi') if re.match('\dpi', x) else x))
        for i in range(NUM_CV):
            header.append(CV(name=firstline[i + 2],
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


def load_colvar(inp_file, cv_names):
    """
    Function to load the colvar colums for the given CVS
    Parameters
    ----------
    inp_file : str
        path to ``inp_file`` (COLVAR file)
    cv_names : list
        List of CV names to load

    Returns
    -------
    colvar : numpy.ndarray
        COLVAR for the given `cv_names`.

    """
    with open(inp_file, 'r') as fp:
        firstline = fp.readline().split()[2:]
    usecols = [0] + [firstline.index(name) for name in cv_names]
    return np.genfromtxt(inp_file, usecols=usecols)


def plot_1D(header, data, fes_edges=None, ax=None):
    # Type: (list, numpy.ndarray) -> matplotlib.image.AxesImage
    """

    Function to plot a 2D ``fes.dat`` file.

    Parameters
    ----------
    header : list
        header of the ``fes.dat`` file
    data : numpy.ndarray
        fes.dat
    fes_edges : List[numpy.ndarray] or None
        List of edges from the fes. If None they will be generated with the min, max, bins of the header.
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
    else:
        fig = ax.get_figure()

    # create fes_edges if not provided
    if fes_edges is None:
        fes_edges = [np.linspace(header[0].min, header[0].max, header[0].nbin)]

    line, = plt.plot(fes_edges[0], data)

    ax.set_xlabel(header[0].name)
    ax.set_ylabel(r'Energy [$\frac{KJ}{mol}$]')

    ax.set_xlim([header[0].min, header[0].max])
    ax.set_ylim([data.min(), data.max()])

    return line


def plot_2D(header, data, ax=None):
    # Type: (list, numpy.ndarray) -> matplotlib.image.AxesImage
    """

    Function to plot a 2D ``fes.dat`` file.

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
    else:
        fig = ax.get_figure()

    extent = list(chain.from_iterable(map(lambda x: [x.min, x.max], header)))

    img = ax.imshow(data, origin='lower',
                    extent=extent,
                    aspect=float(extent[1] - extent[0]) / (extent[3] - extent[2])
                    )

    ax.set_xlabel(header[0].name)
    ax.set_ylabel(header[1].name)

    ax.set_xlim([header[0].min, header[0].max])
    ax.set_ylim([header[1].min, header[1].max])

    cbar = plt.colorbar(img, ax=ax)
    cbar.set_label(r'Energy [$\frac{KJ}{mol}$]')

    return img, cbar


def update(t):
    """
    Update function to move the current point in the graph.

    Parameters
    ----------
    t : int
        frame
    """
    find = colvar_data[colvar_data[:, 0] == t]
    if find.size == 0:
        return
    if colvar_data.shape[1] == 2:
        time, x = find[0]
        y = sorted(zip(abs(fes_edges[0] - x), fes_data),
                   key=lambda x: x[0])[0][1]
    else:
        time, x, y = find[0]
    point.set_data([x], [y])
    point.get_figure().canvas.draw()


def register():
    """
    Function to register VMD function `connect_plot`
    which calls the python function `update($frame)` from within VMD
    """
    evaltcl("""
    proc connect_plot {name index op} {
        # name == vmd_frame
        # index == molecule id of the newly changed frame
        # op == w
        set frame [molinfo $index get frame]
        gopython -command "update($frame)" 
        return
    }
    """)


def start_trace():
    """
    Function to connect VMD function `connect_plot` with VMD variable change of `vmd_frame`
    """
    global fig
    molid = evaltcl("molinfo top")
    evaltcl("trace variable vmd_frame({}) w connect_plot".format(molid))
    fig.canvas.mpl_connect('close_event', stop_trace)


def stop_trace():
    """
    Function to disconnect VMD function `connect_plot` with VMD variable change of `vmd_frame`
    """
    molid = evaltcl("molinfo top")
    evaltcl("trace vdelete vmd_frame({}) w connect_plot".format(molid))


def run(fes_file='fes.dat', colvar='COLVAR', connect=True):
    """
    Run routine to plot the free energy surface and connect the current VMD to a point of it.

    Parameters
    ----------
    fes_file : str
        path to fes file. (Default is 'fes.dat'.)
    colvar : str
        path to COLVAR file (Default is 'COLVAR'.)
    connect : bool
        Switch to connect the plot to the current VMD frame.
    """

    assert os.path.exists(fes_file), "FES File: {} not found".format(fes_file)
    assert os.path.exists(colvar), "COLVAR file: {} not found".format(colvar)

    molid = evaltcl("molinfo top")
    numframes = int(evaltcl('molinfo {} get numframes'.format(molid)))
    current_frame = int(evaltcl('molinfo {} get frame'.format(molid)))

    # load files
    global colvar_data, fes_edges, fes_data
    header, fes_edges, fes_data = load(fes_file)
    dimensions = len(fes_edges)
    colvar_data = load_colvar(colvar, [cv.name for cv in header])

    global fig, ax
    if dimensions == 1:
        # create plot
        fig, ax = plt.subplots()
        line = plot_1D(header, fes_data, fes_edges, ax=ax)
    elif dimensions == 2:
        # create plot
        fig, ax = plt.subplots()
        img, cbar = plot_2D(header, fes_data, ax=ax)
    else:
        raise NotImplementedError("Only 1D and 2D free energy landscapes are implemented." + \
                                  "\nFound a {}D landscape!".format(dimensions))

    # draw
    global point
    point, = ax.plot([0], [0], 'r*', markersize=12)

    try:
        fig.show()
    except AttributeError as e:
        print(e)

    update(current_frame)

    if connect:
        print("register connect_plot")
        register()
        print("Connect plot to frames")
        start_trace()
        print("Stop connect_plot with:\n" +
              "trace vdelete vmd_frame({}) w connect_plot".format(molid))

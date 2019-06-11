import pickle
import subprocess
import sys


def create_standalone_plot(fig, fname, backend=None):
    """
    Create a script which can be executed to plot the given figure.

    Pickles the figure and stores it as string in the script.

    Parameter
    ---------
    fig : matplotlib.figure.Figure
        Matplotlib figure to store.
    fname : str
        File name.
    backend : str or None, optional
        Sets the used backend. Default is None.
        Expamle: 'Qt5Agg, TkAgg'

    Examples
    --------
    Normal
    >>> create_standalone_plot(fig, 'pmf')

    Changing the backend
    >>> create_standalone_plot(fig, 'pmf', backend='Qt5Agg')
    """

    def in_ipynb():
        return 'ipykernel' in sys.modules

    pkl_string = pickle.dumps(fig, protocol=0)

    with open(fname, 'w') as fp:
        fp.write('#!/usr/bin/env python{}.{} \n'.format(sys.version_info.major, sys.version_info.minor))

        fp.write('import pickle \n')
        if backend is not None:
            fp.write('import matplotlib \n')
            fp.write('matplotlib.use("{}")\n'.format(backend))
        fp.write('import matplotlib.pyplot as plt \n')

        if sys.version_info.major < 3:
            import base64
            fp.write("import base64 \n")
            fp.write("pkl_string = b'''{}''' \n".format(base64.b64encode(pkl_string)))
            fp.write('fig = pickle.loads( base64.b64decode(pkl_string) ) \n')

        else:
            fp.write("pkl_string = {} \n".format(pkl_string))
            fp.write('fig = pickle.loads(pkl_string) \n')

        if in_ipynb():
            fp.write('fig._original_dpi = {} \n'.format(fig.get_dpi()))
            fp.write('dummy = plt.figure(figsize={}, dpi={}) \n'.format(
                tuple(fig.get_size_inches()), fig.get_dpi()))
            fp.write('new_manager = dummy.canvas.manager \n')
            fp.write('new_manager.canvas.figure = fig \n')
            fp.write('fig.set_canvas(new_manager.canvas) \n')

        fp.write('plt.show() \n')
    subprocess.Popen("chmod +x {}".format(fname), shell=True)

    print("Created : \033[0;31m{}\033[0m".format(fname))
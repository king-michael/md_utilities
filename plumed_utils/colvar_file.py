# coding=utf-8
"""
Functions around COLVAR files
"""

import os
import numpy as np


def combine_colvar_files(root, colvar_files):
    """
    Combine multiple COLVAR files

    Parameters
    ----------
    root : str
        path where the colvar_files are.
    covlar_files : List[str]
        List of COLVAR files

    Returns
    -------
    new_header : list
        Combined header
    new_data : numpy.ndarray
        Combined data

    Raises
    ------
    AssertionError
        `'Could not find column time'`
         If the COLVAR file does not contain a time column.
    AssertionError
        `'Time steps does not fit'`
        If the COLVAR files does not contain the same time steps.

    Examples
    --------
    >>> root = "path/to/files"
    >>> colvar_files=['COLVAR.distances.dat', 'COLVAR.angles.dat', 'COLVAR.dihedrals.dat']
    >>> header, data = combine_colvar_files(root, colvar_files)
    """

    for i, cfile in enumerate(colvar_files):
        fpath = os.path.join(root, cfile)

        # get header
        with open(fpath, 'r') as fp:
            header = fp.readline().split()[2:]  # get the header of the COLVAR file
        # get the column for time
        i_time = header.index('time')
        assert i_time >= 0, 'Could not find column time'

        # read in data
        data = np.genfromtxt(fpath)

        # get the time column
        time = data[:, i_time]
        header.pop(i_time)
        if i == 0:
            new_data = [np.atleast_2d(time).T]
            new_header = ['time']
        new_header.extend(header)
        np.testing.assert_array_equal(time, new_data[0][:, 0], err_msg='Time steps does not fit')

        new_data.append(np.delete(data, i_time, 1))

    new_data = np.hstack(new_data)

    return new_header, new_data

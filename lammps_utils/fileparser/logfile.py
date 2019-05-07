# coding=utf-8
import re
import numpy as np


def read_logfile(fname, combine=True, as_recarray=True, try_corrupt_file=True):
    """
    Function to read a LAMMPS log file.

    Notes
    -----
    * ``thermo_style multi`` is not supported
    * Every column is read as ``np.float64``.  Maximum integer without rounding:
        * ``np.float32`` : ``16777216``
        * ``np.float64`` : ``9007199254740992``

    Parameters
    ----------
    fname : str
        path to the LAMMPS log file.
    combine : bool, optional
        Switch if different runs in a logfile should be combined into one or not
    as_recarray : bool, optional
        Switch if the result should be returned as ``numpy.recarray``.
        The column names will be used as field names.
    try_corrupt_file : bool, optional
        Try to read a corrupt file where the run didn't finish.

    Returns
    -------
    header : List[str] or List[List[str]]
        List of columns. If ``combine`` is ``False`` a list of lists will be returned.
    data : np.ndarray or np.recarray or List[np.ndarray] or List[np.recarray]
        Extracted data. If ``combined`` is ``False`` a list of numpy arrays will be returned.
        If ``as_recarray`` is ``False``, the data is stored in a ``np.ndarray``
        otherwise in a ``np.recarray`` with the column names from the header as field names.

    Raises
    ------
    AssertionError
        ``'No thermodynamic output found in file : {}'``
        If no thermodynamic output is found.
    AssertionError
        ``'Different thermodynamic output can not combine it'``
        If ``combine=True`` but the different runs have different thermodynamic properties.

    Examples
    --------

    Default

    >>> header, data = read_logfile('log.lammps')
    >>> # test if data is a np.recarray
    >>> type(data) == np.recarray
    True
    >>> # test if the array has the same names as columns in the header
    >>> header == list(data.dtype.names)
    True
    Return the values as numpy array

    >>> header, recarray = read_logfile('log.lammps', as_recarray=False)
    >>> # test if data is a np.ndarray
    >>> type(data) == np.ndarray
    True

    Read a logfile with multiple runs in it and do not combine them.

    >>> header, data = read_logfile('log.lammps', combine=False)
    >>> # test if data is a list and the entry is a np.recarray
    >>> type(data) == list and type(data[0]) == np.recarray
    True
    >>> # test if all arrays have the same names as columns in the header
    >>> all(h == list(d.dtype.names) for h, d in zip(header, data))
    True

    """

    pattern_data = re.compile('(^Step.*$)\n((?s:.*))\nLoop', re.MULTILINE)
    pattern_data_corrupt = re.compile('(^Step.*$)\n((?s:.*))', re.MULTILINE)

    dtype = np.float64  # used dtype for the columns

    list_header = []
    list_data = []

    with open(fname, 'r') as fp:
        text = fp.read()

    for subtext in re.finditer(pattern_data, text):
        list_header.append(subtext.group(1).split())
        list_data.append(np.genfromtxt(subtext.group(2).split('\n'), dtype=dtype))

    if try_corrupt_file:
        # get the position to start from
        pos_start = subtext.end() if len(list_header) > 0 else 0
        # extract data
        for subtext in re.finditer(pattern_data_corrupt, text[pos_start:]):
            list_header.append(subtext.group(1).split())
            list_data.append(np.genfromtxt(subtext.group(2).split('\n'), dtype=dtype))
            # TODO : add handling if line is not complete

    assert len(list_header) > 0, 'No thermodynamic output found in file : {}'.format(fname)

    if len(list_header) == 1:
        # If it is only one data block
        header = list_header[0]
        data = list_data[0]

        if as_recarray:  # convert to recarray
            data = np.rec.fromarrays(data.T, dtype=np.dtype([(f, dtype) for f in header]))

    elif combine:
        # If there are serveral data blocks and they should be combined
        assert all([h == list_header[0] for h in list_header]), 'Different thermodynamic output can not combine it'

        header = list_header[0]
        data = np.vstack(list_data)

        if as_recarray:  # convert to recarray
            data = np.rec.fromarrays(data.T, dtype=np.dtype([(f, dtype) for f in header]))

    else:
        # If there are serveral data blocks and they should not be combined
        header = list_header
        data = list_data

        if as_recarray:  # convert to recarray
            data = [np.rec.fromarrays(array.T, dtype=np.dtype([(f, dtype) for f in header]))
                    for header, array in zip(list_header, data)]

    return header, data

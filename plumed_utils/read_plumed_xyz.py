import numpy as np

def read_plumed_xyz(fname, len_element=2):
    """
    Generator to read a PLUMED xyz file.

    Parameters:
    -----------
    fname : str
        path to file

    Yield:
    ------
    data: numpy.ndarray
        array of the data of the corresponding timestep.

    Example:
    --------
    Get all data:
    >>> data = np.array([data for data in read_dump(filename)])
    Get only forces:
    >>> data = np.array([np.transpose([data['fx'], data['fy'], data['fz']])
    >>>                  for data in read_dump(filename)])
    """

    if fname[-3:] == '.gz':
        import gzip
        fp = gzip.open(fname, 'r')
        readline = lambda : fp.readline().decode()
    else:
        fp = open(fname, 'r')
        readline = lambda : fp.readline()

    dtype_str='<U{:d}'.format(len_element)
    line = ' '
    while line:
        line = readline().strip()
        if len(line) == 0:
            continue
        n_atoms = int(line)
        box = tuple(readline().strip().split())
        # read in data
        line_split = readline().split()
        # number of attributes
        n_attr = len(line_split) - 1 - 3
        # initialize data
        element = np.empty((n_atoms, 1), dtype=dtype_str)
        xyz = np.empty((n_atoms, 3), dtype=np.float)
        data = np.empty((n_atoms, n_attr), dtype=np.float)
        # set first row
        element[0], xyz[0], data[0] = (line_split[0],
                                       tuple(line_split[1:4]),
                                       tuple(line_split[4:]))
        # iterate over the file
        for i in range(n_atoms-1):
            line_split = readline().split()
            element[i], xyz[i], data[i] = (line_split[0],
                                           tuple(line_split[1:4]),
                                           tuple(line_split[4:]))
        yield element, xyz, data
    fp.close()
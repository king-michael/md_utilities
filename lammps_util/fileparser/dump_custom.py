#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Parser for the custom dump file.

Author : Michael King <michael.king@uni-konstanz.de>
"""

import numpy as np


def read_dump(fname):
    """
    Generator to read dump files.

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
    >>> data = np.array([data for data in read_dump(filename)])
    """

    if fname[-3:] == '.gz':
        import gzip
        fp = gzip.open(fname, 'r')
    else:
        fp = open(fname, 'r')

    line = ' '
    while line:
        line = fp.readline()

        line_decode = line.decode().strip()
        line_split = line_decode.split()
        if line_decode == 'ITEM: TIMESTEP':
            ts = int(fp.readline().strip())
        elif line_decode == 'ITEM: NUMBER OF ATOMS':
            n_atoms = int(fp.readline().strip())
        elif line_decode.startswith('ITEM: BOX BOUNDS'):
            if len(line_split) == 6:
                triclinic = False
            elif len(line_split) == 9:
                triclinic = True
            else:
                raise NotImplementedError('Not implemented this type of a box:\n{}'.format(line_decode))
        elif line_decode.startswith('ITEM: ATOMS'):
            properties = line_split[2:]
            data = np.zeros(n_atoms,
                            dtype=[(prop, np.uint32 if prop in ['id', 'type'] else np.float)
                                   for prop in properties])
            for i in range(n_atoms):
                data[i] = tuple(fp.readline().decode().split())
            yield data

#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Parser for the custom dump file.

Author : Michael King <michael.king@uni-konstanz.de>

Todo
----
implement dtype
    setting, so we can define which dtype our output should have (decreases memory)
implement n_atoms_constant = True:
    then we can try to build a function that:
        1.) parse the whole script,
        2.) gets number of lines, c
        3.) calculates number timesteps in the file
        4.) initalize the data array
        5.) fast parse the file
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
    Get all data:
    >>> data = np.array([data for data in read_dump(filename)])
                  for data in read_dump(filename)])
    """
    #Get only forces:
    #>>> data = np.array([np.transpose([data['fx'], data['fy'], data['fz']])
    #>>>

    if fname[-3:] == '.gz':
        import gzip
        fp = gzip.open(fname, 'r')
        readline = lambda: fp.readline().decode()
    else:
        fp = open(fname, 'r')
        readline = lambda: fp.readline()

    line = ' '
    while line:
        line = readline().strip()

        if line == 'ITEM: TIMESTEP':  # "ITEM: TIMESTEP\n"
            ts = int(readline().strip())
            line = fp.readline()  # "ITEM: NUMBER OF ATOMS\n"
            n_atoms = int(readline().strip())

            # orthogonal: "ITEM: BOX BOUNDS %s\n",boundstr
            # triclinic: "ITEM: BOX BOUNDS xy xz yz %s\n",boundstr
            line_decode = readline().strip()
            line_split = line_decode.split()

            if len(line_split) == 6: # orthogonal box
                # "%-1.16e %-1.16e %-1.16e\n",boxxlo,boxxhi,boxxy
                boxxlo, boxxhi, boxxy = readline().strip().split()
                boxylo, boxyhi, boxxz = readline().strip().split()
                boxzlo, boxzhi, boxyz = readline().strip().split()

            if len(line_split) == 9: # triclinic box
                # "%-1.16e %-1.16e\n",boxxlo,boxxhi
                boxxlo, boxxhi = readline().strip().split()
                boxylo, boxyhi = readline().strip().split()
                boxzlo, boxzhi = readline().strip().split()

            line_decode = readline().strip()  # "ITEM: ATOMS %s\n",columns
            line_split = line_decode.split()
            properties = line_split[2:]
            # data = np.zeros(n_atoms,
            #                 dtype=[(prop, np.uint32 if prop in ['id', 'type'] else np.float)
            #                        for prop in properties])
            data = np.zeros((n_atoms, len(properties)), dtype=np.float)

            for i in range(n_atoms):
                data[i] = tuple(readline().split())
            yield data

    fp.close()

#!/usr/bin/env python
# coding=utf-8
"""
Various Utilities around LAMMPS
"""

import os
import re


def extract_unit_constants_from_src(path_to_lammps_dir):
    """
    Function to extract the constants used in LAMMPS formulas

    Parameters
    ----------
    path_to_lammps_dir : str
        Path to the LAMMPS root directory.

    Returns
    -------
    constants : dict
        Dictionary with the constants for every unit.
    """
    file_path = os.path.join(path_to_lammps_dir, 'src', 'update.cpp')

    with open(file_path) as fp:
        text = fp.read()

    pattern_styles = re.compile(r'if \(strcmp\(style,"(\w+)"\) == 0\) {[\w\d ;\n->=]*}', re.MULTILINE)
    pattern_constants = re.compile(r'force->(\w+) = (.*);')

    constants = dict(
        (
            m.group(1),
            dict(
                (c.group(1), eval(c.group(2)))
                for c in re.finditer(pattern_constants, text[slice(*m.span())]))
        ) for m in re.finditer(pattern_styles, text)
    )

    return constants

#!/usr/bin/env python
"""
Functions to get thermodynamic properties from a logfile.
"""


from __future__ import absolute_import

import numpy as np
from .logfile import read_logfile, get_dof, get_n_atoms


def get_cv(logfile, temperature=None):
    """
    Extract the isochoric heat capacity from the fluctuation of the total energy in a LAMMPS logfile.

    Needs `units real`

    Parameter
    ---------
    logfile : str
        path to the LAMMPS logfile
    temperature : float or None, optional
        Temperature to use. If `None` the mean temperature of the simulation is used.
        Default is `None`.
    Returns
    -------
    cv : float
        Isochoric heat capacity [J/(mol * K)]
    """
    header, data = read_logfile(logfile)

    # R = 8.314462618 /1000 / 4.184 # J/K/mol # kJ/K/mol -> kcal/K/mol
    boltz = 0.0019872067  # kcal/K/mol value from lammps

    temp = data['Temp']

    if 'KinEng' in header:
        ke = data['KinEng']
    else:
        n_dof = get_dof(logfile)
        ke = n_dof / 2 * boltz * data['Temp']

    if 'TotEng' in header:
        etotal = data['TotEng']
    else:
        etotal = data['PotEng'] + ke

    if temperature is None:
        temperature = np.mean(temp)

    # np.var(etotal*1000) / (boltz * 1000 * np.power(temperature,2)) * 4.184 / 1000
    # np.var(etotal*1000 * 4.184) / (8.314462618 * np.power(temperature,2)) / 1000
    n_atoms = get_n_atoms(logfile)
    cv = np.var(etotal) / (boltz * np.power(temperature, 2)) * 4.184

    return cv


def get_cp(logfile, temperature=None, pressure='mean', atoms_per_molecule=1):
    """
    Extract the isobaric heat capacity from the fluctuation of the total energy in a LAMMPS logfile.

    Needs `units real`

    Parameter
    ---------
    logfile : str
        path to the LAMMPS logfile
    temperature : float or None, optional
        Temperature to use. If `None` the mean temperature of the simulation is used.
        Default is `None`.
    Returns
    -------
    cp : float
        Isobaric heat capacity [J/(mol * K)]
    """
    header, data = read_logfile(logfile)

    # R = 8.314462618 /1000 / 4.184 # J/K/mol # kJ/K/mol -> kcal/K/mol
    boltz = 0.0019872067  # kcal/K/mol value from lammps

    temp = data['Temp']

    if 'KinEng' in header:
        ke = data['KinEng']
    else:
        n_dof = get_dof(logfile)
        ke = n_dof / 2 * boltz * data['Temp']

    if 'TotEng' in header:
        etotal = data['TotEng']
    else:
        etotal = data['PotEng'] + ke

    if temperature is None:
        temperature = np.mean(temp)

    if type(pressure) is str:
        if pressure == 'mean':
            pressure = np.mean(data['Press'])
        else:
            pressure = data['Press']

    if 'Enthalpy' in header:
        enthalpy = data['Enthalpy']
    else:
        conv_vp2e = 1.0 / 68568.415
        enthalpy = etotal + pressure * data['Volume'] * conv_vp2e

    n_atoms = get_n_atoms(logfile)
    # enthalpy *= 1000
    R = 8.314462618
    cp = np.var(enthalpy) / (n_atoms / atoms_per_molecule) / (boltz * np.power(temperature, 2)) * 1000 * 4.184

    return cp
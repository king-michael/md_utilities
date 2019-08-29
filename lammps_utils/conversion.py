# coding=utf-8
"""
Different functions to convert simulation properties

Functions
---------
* ``convert_box`` : convert different box formats

"""
import numpy as np


def convert_box(box,
                return_box_vectors=False,
                return_box_lammps=False,
                triclinic=False,
                precision=15):
    """
    Function to convert box formats

    Parameters
    ----------
    box : list, List[list], np.ndarray
        Simulation box can have the format:
        * normal : ``[a, b, c, alpha, beta, gamma]``,
        * box_vectors : ``[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]``,
        * LAMMPS orthogonal : ``[[xlo, xhi], [ylo, yhi], [zlo, zhi]]``,
        * LAMMPS triclinic : ``[[xlo, xhi], [ylo, yhi], [zlo, zhi], [xy, xz, yz]]``
    return_box_vectors : bool, optional
        returns box vectors in format ``[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]``
    return_box_lammps : bool, optional
        returns box in lammps format.
        If all angles are 90 degree it returns a orthogonal box with ``[[xlo, xhi], [ylo, yhi], [zlo, zhi]]``.
        Otherwise or if ``triclinic=True`` it returns a box with ``[[xlo, xhi], [ylo, yhi], [zlo, zhi], [xy, xz, yz]]``.
    triclinic : bool, optional
        together with ``return_box_lammps`` returns triclinic box even when all angles are 90 degree
    precision : int, optional
        precision for rounding

    Returns
    -------
    box : list or List[list]
        normal : ``[a, b, c, alpha, beta, gamma]``
        box_vectors : ``[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]``
        box_lammps(orthogonal) : ``[[xlo, xhi], [ylo, yhi], [zlo, zhi]]``
        box_lammps(triclinic) : ``[[xlo, xhi], [ylo, yhi], [zlo, zhi], [xy, xz, yz]]``

    """
    calculate_box = False
    calculate_box_vectors = False
    calculate_box_bounds = False
    if len(box) == 6:
        a, b, c, alpha, beta, gamma = box

        calculate_box_vectors = True
        calculate_box_bounds = True

    elif len(box) == 3:
        if 3 == len(box[0]) == len(box[1]) == len(box[2]):
            box_vectors = box

            (lx, _, _), (xy, ly, _), (xz, yz, lz) = box_vectors

            calculate_box = True
            calculate_box_bounds = True

        elif 2 == len(box[0]) == len(box[1]) == len(box[2]):
            (xlo, xhi), (ylo, yhi), (zlo, zhi) = box
            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo
            xy = xz = yz = 0.0
            calculate_box = True
        else:
            raise UserWarning("Box could not be read.\nbox={}".format(box))

    elif len(box) == 4:
        if 2 == len(box[0]) == len(box[1]) == len(box[2]) and len(box[3]) == 3:
            (xlo, xhi), (ylo, yhi), (zlo, zhi), (xy, xz, yz) = box
            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo
            calculate_box = True

        else:
            raise UserWarning("Box could not be read.\nbox={}".format(box))

    else:
        raise UserWarning("Box could not be read.\nbox={}".format(box))

    if calculate_box:
        a = lx
        b = np.sqrt(ly ** 2 + xy ** 2)
        c = np.sqrt(lz ** 2 + xz ** 2 + yz ** 2)
        alpha = np.round(np.rad2deg(np.arccos((xy * xz + ly * yz) / (b * c))), precision)
        beta = np.round(np.rad2deg(np.arccos(xz / c)), precision)
        gamma = np.round(np.rad2deg(np.arccos(xy / b)), precision)

    if calculate_box_vectors:
        lx = a
        xy = b * np.round(np.cos(np.deg2rad(gamma)), precision)
        xz = c * np.round(np.cos(np.deg2rad(beta)), precision)
        ly = np.sqrt(b ** 2 - xy ** 2)
        yz = (b * c * np.round(np.cos(np.deg2rad(alpha)), precision) - xy * xz) / ly
        lz = np.sqrt(c ** 2 - xy ** 2 - yz ** 2)

    if calculate_box_bounds:
        xlo, xhi = 0.0, lx
        ylo, yhi = 0.0, ly
        zlo, zhi = 0.0, lz

    if return_box_lammps:
        if triclinic or not xy == xz == yz == 0:
            return [[xlo, xhi], [ylo, yhi], [zlo, zhi], [xy, xz, yz]]

        return [[xlo, xhi], [ylo, yhi], [zlo, zhi]]

    if return_box_vectors:
        return [[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]

    return [a, b, c, alpha, beta, gamma]
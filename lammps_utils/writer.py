#!/usr/bin/env python
# coding=utf-8
import numpy as np
import sys
import datetime


class LAMMPSDATA:
    def __init__(self,
                 positions,
                 box,
                 types=None,
                 ids=None,
                 masses=None,
                 first_line=None,
                 ):
        """
        Writer for LAMMPS data files
        Parameters
        ----------
        positions
        box
        types
        ids
        masses
        first_line
        """
        # @TODO: add fields (type_names, names, bonds, angles, dihedrals, impropers, etc)
        self.positions = positions
        self.box = box
        self.types = types
        self.ids = ids
        self.masses = masses

        self.first_line = first_line
        self._check_integrity()

    def _check_integrity(self):
        """checks if all parameter are ok and does some conversions"""
        if self.first_line is not None:
            self.first_line = self.first_line.rstrip()
            assert self.first_line.find('\n') == -1, \
                "Only one line allowed for 'first_line'\nfirst_line='{}'".format(self.first_line)

        self.positions = np.asarray(self.positions)
        self.n_atoms = self.positions.shape[0]

        if self.types is None:
            self.types = np.ones(self.n_atoms, dtype=int)
        else:
            self.types = np.asarray(self.types)
            assert self.types.size == self.n_atoms, \
                "Number of atoms does not match number of types"

        self.types_uniq = sorted(set(self.types))
        self.n_atom_types = len(self.types_uniq)
        self.atype_names_uniq = [str(t) for t in self.types_uniq]  # DUMMY
        # @TODO : types names to int

        if self.ids is None:
            self.ids = np.arange(self.n_atoms, dtype=int) + 1
        else:
            self.ids = np.asarray(ids)
            assert self.ids.size == self.n_atoms, \
                "Number of atoms does not match number of ids"

        if self.masses is None:
            self.masses = np.ones(self.n_atom_types)
        else:
            self.masses = np.asarray(masses)
            assert self.masses.size == self.n_atom_types, \
                "Number of atoms does not match number of masses"

        self.n_bonds = 0  # DUMMY
        self.n_bond_types = 0  # DUMMY
        self.n_angles = 0  # DUMMY
        self.n_angle_types = 0  # DUMMY
        self.n_dihedrals = 0  # DUMMY
        self.n_dihedral_types = 0  # DUMMY
        self.n_impropers = 0  # DUMMY
        self.n_improper_types = 0  # DUMMY

        self.charges = np.zeros(self.n_atoms) # DUMMY
        self.names = ['X' for _ in range(self.n_atoms)]  # DUMMY
        self.resnames = ['X' for _ in range(self.n_atoms)]  # DUMMY
        self.resid = np.zeros(self.n_atoms, dtype=int)  # DUMMY

    def convert_box(self,
                    box,
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

    def write(self, fname=None):
        """
        Writes LAMMPS data file.

        Parameters
        ----------
        fname : str or None, optional
            Name of the datefile. If ``None`` it will be writen to the ``sys.stdout``.
            Default is ``None``.
        """
        fp = sys.stdout if fname is None else open(fname, 'w')
        if self.first_line is None:
            fp.write(('LAMMPS data file. atom_style full generated by Python on '
                      '{:%Y-%m-%d %H:%M} UTC\n').format(datetime.datetime.now(datetime.timezone.utc)))

        else:
            fp.write('{first_line}\n'.format(self.first_line))
        fp.write((
                     ' {n_atoms:d} atoms\n'
                     ' {n_bonds:d} bonds\n'
                     ' {n_angles:d} angles\n'
                     ' {n_dihedrals:d} dihedrals\n'
                     ' {n_impropers:d} impropers\n'
                     ' {n_atom_types:d} atom types\n'
                     ' {n_bond_types:d} bond types\n'
                     ' {n_angle_types:d} angle types\n'
                     ' {n_dihedral_types:d} dihedral types\n'
                     ' {n_improper_types:d} improper types\n'
                 ).format(**self.__dict__))
        # BOX
        box_lammps = self.convert_box(self.box, return_box_lammps=True)
        fp.write((' {:g} {:g} xlo xhi\n'
                  ' {:g} {:g} ylo yhi\n'
                  ' {:g} {:g} zlo zhi\n').format(*box_lammps[0], *box_lammps[1], *box_lammps[2]))
        if len(self.box) == 4:
            fp.write(' {:g} {:g} {:g} xy xz yz\n'.format())

        # Masses
        fp.write('\n Masses\n\n')
        fp.write('\n'.join([' {} {} # {}'.format(t, m, n)
                            for m, t, n in zip(self.masses, self.types_uniq, self.atype_names_uniq)]))
        fp.write('\n')

        # Atoms
        fp.write('\n Atoms\n\n')
        # TODO : implement other line
        line_format = '{id:d} {resid:d} {type:d} {charge:f} {x:f} {y:f} {z:f} # {name:%ds} {resname:%ds}\n'
        line_format = line_format % (max(len(s) for s in self.names), max(len(s) for s in self.resnames))  # max len
        for i, t, r, c, x, y, z, n, rn in zip(
                self.ids, self.types, self.resid, self.charges,
                self.positions[:, 0], self.positions[:, 1], self.positions[:, 2],
                self.names, self.resnames):
            fp.write(line_format.format(id=i, type=t, resid=r, charge=c,
                                        x=x, y=y, z=z, name=n, resname=rn))

        # @TODO : Bonds
        # @TODO : Angles
        # @TODO : Diheadrals
        # @TODO : Impropers

        if fname is not None:
            fp.close()
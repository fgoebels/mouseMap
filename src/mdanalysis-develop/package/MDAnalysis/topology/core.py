# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Common functions for topology building --- :mod:`MDAnalysis.topology.core`
==========================================================================

The various topology parsers make use of functions and classes in this
module. They are mostly of use to developers.

.. SeeAlso:: :mod:`MDAnalysis.topology.tables` for some hard-coded atom
   information that is used by functions such as :func:`guess_atom_type` and
   :func:`guess_atom_mass`.

"""

from __future__ import print_function
# Global imports
import os.path
import numpy as np
from collections import defaultdict
from itertools import izip

# Local imports
from . import tables
from ..lib import distances
from ..lib.util import cached
from ..lib import util
from ..core import AtomGroup


def build_segments(atoms):
    """Create all :class:`~MDAnalysis.core.AtomGroup.Segment` instancess from
    a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The function also builds the :class:`~MDAnalysis.core.AtomGroup.Residue`
    instances by tracking residue numbers.

    Updating segments also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue and segment an atom belongs.

    :Returns: List of :class:`~MDAnalysis.core.AtomGroup.Segment`
              instances

    .. versionchanged:: 0.9.0
       Now allows resids in a given :class:`Segment` to be given in non sequential order.
    .. versionchanged:: 0.11.0
       Returns a list instead of a dict for consistency with the outputs of
       :func:`build_residues`
    """
    struc = list()
    resatomlist = defaultdict(list)
    curr_segname = atoms[0].segid

    for a in atoms:
        if a.segid == curr_segname:  # if still in same Segment
            resatomlist[a.resid].append(a)
        else:
            # We've come to a new segment
            # Build the Segment we just left
            residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                        for k, ats in resatomlist.iteritems()]
            struc.append(AtomGroup.Segment(curr_segname, residues))

            # Reset things and start again
            resatomlist = defaultdict(list)
            resatomlist[a.resid].append(a)
            curr_segname = a.segid

    # Add the last segment
    residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                for k, ats in resatomlist.iteritems()]
    struc.append(AtomGroup.Segment(curr_segname, residues))
    return struc


def build_residues(atoms):
    """Create a list :class:`~MDAnalysis.core.AtomGroup.Residue` instances from
    a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Updating residues also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue an atom belongs.

    :Returns: List of :class:`~MDAnalysis.core.AtomGroup.Residue` instances

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Now allows resids to be given in non sequential order
    """
    resatomlist = defaultdict(list)
    for a in atoms:
        resatomlist[a.resid].append(a)

    residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                for k, ats in resatomlist.iteritems()]

    return residues


def get_parser_for(filename, permissive=False, format=None):
    """Return the appropriate topology parser for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.

    :Raises:
      *ValueError*
        If no appropriate parser could be found.
    """
    from . import _topology_parsers

    if format is None:
        format = util.guess_format(filename)
    format = format.upper()
    if format == 'PDB' and permissive:
        return _topology_parsers['Permissive_PDB']

    try:
        return _topology_parsers[format]
    except KeyError:
        raise ValueError(
            "Cannot autodetect topology type for file '{0}' "
            "(file extension could not be parsed).\n"
            "           You can use 'Universe(topology, ..., topology_format=FORMAT)' "
            "to explicitly specify the format and\n"
            "           override automatic detection. Known FORMATs are:\n"
            "           {1}\n"
            "           See http://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats\n"
            "           For missing formats, raise an issue at "
            "http://issues.mdanalysis.org".format(filename, _topology_parsers.keys()))


# following guess_* used by PDB parser
def guess_atom_type(atomname):
    """Guess atom type from the name.

    At the moment, this function simply returns the element, as
    guessed by :func:`guess_atom_element`.

    .. SeeAlso:: :func:`guess_atom_element` and
                 :mod:`MDAnalysis.topology.tables`
    """
    return guess_atom_element(atomname)


def guess_atom_element(atomname):
    """Guess the element of the atom from the name.

    Looks in dict to see if element is found, otherwise it uses the first
    character in the atomname. The table comes from CHARMM and AMBER atom
    types, where the first character is not sufficient to determine the atom
    type. Some GROMOS ions have also been added.

    .. Warning: The translation table is incomplete. This will probably result
                in some mistakes, but it still better than nothing!

    .. SeeAlso:: :func:`guess_atom_type` and
                 :mod:`MDAnalysis.topology.tables` (where the data are stored)
    """
    try:
        return tables.atomelements[atomname]
    except KeyError:
        if atomname[0].isdigit():
            # catch 1HH etc
            return atomname[1]
        return atomname[0]


def guess_bonds(atoms, coords, **kwargs):
    """Guess if bonds exist between two atoms based on their distance.

    Bond between two atoms is created, if the two atoms are within

    .. math::

          d < f * (R_1 + R_2)

    of each other, where :math:`R_1` and :math:`R_2` are the VdW radii
    of the atoms and :math:`f` is an ad-hoc *fudge_factor*. This is
    the `same algorithm that VMD uses`_.

    :Keywords:

      *fudge_factor*
        The factor by which atoms must overlap eachother to be considered a
        bond.  Larger values will increase the number of bonds found. [0.72]

      *vdwradii*
        To supply custom vdwradii for atoms in the algorithm. Must be a dict
        of format {type:radii}. The default table of van der Waals radii is
        hard-coded as :data:`MDAnalysis.topology.tables.vdwradii`.  Any user
        defined vdwradii passed as an argument will supercede the table
        values. [``None``]

      *lower_bound*
        The minimum bond length. All bonds found shorter than this length will
        be ignored. This is useful for parsing PDB with altloc records where
        atoms with altloc A and B maybe very close together and there should be
        no chemical bond between them. [0.1]

      *box*
        Bonds are found using a distance search, if unit cell information is
        given, periodic boundary conditions will be considered in the distance
        search. [``None``]

    :Returns:
       List of tuples suitable for use in Universe topology building.

    .. warning::
       No check is done after the bonds are guessed to see if Lewis
       structure is correct. This is wrong and will burn somebody.

    :Raises:
       ValueError if inputs are malformed or *vdwradii* data is missing.

    .. _`same algorithm that VMD uses`:
       http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html

    .. versionadded:: 0.7.7
    .. versionchanged:: 0.9.0
       Updated method internally to use more :mod:`numpy`, should work
       faster.  Should also use less memory, previously scaled as
       :math:`O(n^2)`.  *vdwradii* argument now augments table list
       rather than replacing entirely.
    """
    # why not just use atom.positions?
    if len(atoms) != len(coords):
        raise ValueError("'atoms' and 'coord' must be the same length")

    fudge_factor = kwargs.get('fudge_factor', 0.72)

    vdwradii = tables.vdwradii.copy()  # so I don't permanently change it
    user_vdwradii = kwargs.get('vdwradii', None)
    if user_vdwradii:  # this should make algo use their values over defaults
        vdwradii.update(user_vdwradii)

    try:
        atomtypes = set(atoms.types)
    except AttributeError:  # sometimes atoms is just list of atoms not AG
        atomtypes = set([a.type for a in atoms])
    # check that all types have a defined vdw
    if not all([val in vdwradii for val in atomtypes]):
        raise ValueError(("vdw radii for types: " +
                          ", ".join([t for t in atomtypes if
                                     not t in vdwradii]) +
                          ". These can be defined manually using the" +
                          " keyword 'vdwradii'"))

    lower_bound = kwargs.get('lower_bound', 0.1)

    box = kwargs.get('box', None)

    # to speed up checking, calculate what the largest possible bond
    # atom that would warrant attention.
    # then use this to quickly mask distance results later
    max_vdw = max([vdwradii[t] for t in atomtypes])

    bonds = []

    for i, atom in enumerate(atoms[:-1]):
        vdw_i = vdwradii[atom.type]
        max_d = (vdw_i + max_vdw) * fudge_factor

        # using self_distance_array scales O(n^2)
        # 20,000 atoms = 1.6 Gb memory
        dist = distances.distance_array(coords[i][None, :], coords[i + 1:],
                                        box=box)[0]
        idx = np.where((dist > lower_bound) & (dist <= max_d))[0]

        for a in idx:
            atom_j = atoms[i + 1 + a]

            if dist[a] < (vdw_i + vdwradii[atom_j.type]) * fudge_factor:
                # because of method used, same bond won't be seen twice,
                # so don't need to worry about duplicates
                bonds.append((atom.index, atom_j.index))

    return tuple(bonds)


def guess_angles(bonds):
    """Given a list of Bonds, find all angles that exist between atoms.

    Works by assuming that if atoms 1 & 2 are bonded, and 2 & 3 are bonded,
    then (1,2,3) must be an angle.

    :Returns:
      List of tuples defining the angles.
      Suitable for use in u._topology

    .. seeAlso:: :meth:`guess_bonds`

    .. versionadded 0.9.0
    """
    angles_found = set()

    for b in bonds:
        for atom in b:
            other_a = b.partner(atom)  # who's my friend currently in Bond
            for other_b in atom.bonds:
                if other_b != b:  # if not the same bond I start as
                    third_a = other_b.partner(atom)
                    desc = tuple([other_a.index, atom.index, third_a.index])
                    if desc[0] > desc[-1]:  # first index always less than last
                        desc = desc[::-1]
                    angles_found.add(desc)

    return tuple(angles_found)


def guess_dihedrals(angles):
    """Given a list of Angles, find all dihedrals that exist between atoms.

    Works by assuming that if (1,2,3) is an angle, and 3 & 4 are bonded,
    then (1,2,3,4) must be a dihedral.

    :Returns:
      List of tuples defining the dihedrals.
      Suitable for use in u._topology

    .. versionadded 0.9.0
    """
    dihedrals_found = set()

    for b in angles:
        a_tup = tuple([a.index for a in b])  # angle as tuple of numbers
        # if searching with b[0], want tuple of (b[2], b[1], b[0], +new)
        # search the first and last atom of each angle
        for atom, prefix in zip([b.atoms[0], b.atoms[-1]],
                                [a_tup[::-1], a_tup]):
            for other_b in atom.bonds:
                if not other_b.partner(atom) in b:
                    third_a = other_b.partner(atom)
                    desc = prefix + (third_a.index,)
                    if desc[0] > desc[-1]:
                        desc = desc[::-1]
                    dihedrals_found.add(desc)

    return tuple(dihedrals_found)


def guess_improper_dihedrals(angles):
    """Given a list of Angles, find all improper dihedrals that exist between
    atoms.

    Works by assuming that if (1,2,3) is an angle, and 2 & 4 are bonded,
    then (2, 1, 3, 4) must be an improper dihedral.
    ie the improper dihedral is the angle between the planes formed by
    (1, 2, 3) and (1, 3, 4)

    :Returns:
      List of tuples defining the improper dihedrals.
      Suitable for use in u._topology

    .. versionadded 0.9.0
    """
    dihedrals_found = set()

    for b in angles:
        atom = b[1]  # select middle atom in angle
        # start of improper tuple
        a_tup = tuple([b[a].index for a in [1, 2, 0]])
        # if searching with b[1], want tuple of (b[1], b[2], b[0], +new)
        # search the first and last atom of each angle
        for other_b in atom.bonds:
            other_atom = other_b.partner(atom)
            # if this atom isn't in the angle I started with
            if not other_atom in b:
                desc = a_tup + (other_atom.index,)
                if desc[0] > desc[-1]:
                    desc = desc[::-1]
                dihedrals_found.add(desc)

    return tuple(dihedrals_found)


def get_atom_mass(element):
    """Return the atomic mass in u for *element*.

    Masses are looked up in :data:`MDAnalysis.topology.tables.masses`.

    .. Warning:: Unknown masses are set to 0.00
    """
    try:
        return tables.masses[element]
    except KeyError:
        return 0.000


def guess_atom_mass(atomname):
    """Guess a mass based on the atom name.

    :func:`guess_atom_element` is used to determine the kind of atom.

    .. warning:: Anything not recognized is simply set to 0; if you rely on the
                 masses you might want to double check.
    """
    return get_atom_mass(guess_atom_element(atomname))


def guess_atom_charge(atomname):
    """Guess atom charge from the name.

    .. Warning:: Not implemented; simply returns 0.
    """
    # TODO: do something slightly smarter, at least use name/element
    return 0.0
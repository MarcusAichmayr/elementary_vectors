r"""
Chirotopes constructed from (co)circuits.

EXAMPLES::

    sage: from sign_vectors import *
    sage: from sign_vectors.chirotope import *
    sage: c = ChirotopeFromCircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
    sage: c.entry([0, 1])
    +
    sage: c.entry([0, 2])
    0
    sage: c.entry([1, 2])
    0
    sage: c = ChirotopeFromCircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
    sage: c.entries()
    [+, 0, 0, 0, 0, 0]

::

    sage: c = ChirotopeFromCocircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
    sage: c.entries()
    [0, 0, 0, 0, 0, +]

::

    sage: c = ChirotopeFromCircuits([sign_vector("++0--"), sign_vector("+0--0"), sign_vector("-0++0"), sign_vector("0--0+"), sign_vector("--0++"), sign_vector("0++0-")], 3, 5)
    sage: c.entries()
    [+, -, +, 0, -, +, +, 0, -, +]

::

    sage: c = ChirotopeFromCocircuits([sign_vector("00++"), sign_vector("0++0"), sign_vector("0+0-")], 2, 4)
    sage: c.entries()
    [0, 0, 0, +, +, +]

TESTS:

When iterating over ``{7, 8}``, we obtain ``8`` first.
Therefore, we need to sort indices::

    sage: M = matrix.ones(1, 9)
    sage: om = OrientedMatroid(M)
    sage: ChirotopeFromCircuits(om.circuits(), 1, 9).entries()
    [+, +, +, +, +, +, +, +, +]
"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from enum import IntEnum
from typing import Iterator

from sage.combinat.combination import Combinations

from sign_vectors import sign_symbolic, SignVector


class Sign(IntEnum):
    r"""
    Auxiliary class for chirotopes.

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import Sign
        sage: Sign(1)
        +
        sage: Sign(-1)
        -
        sage: Sign(0)
        0
        sage: Sign(5)
        +
        sage: Sign(5).value
        1
        sage: -Sign(5)
        -
        sage: Sign("+")
        +
        sage: Sign("-")
        -
        sage: Sign("0")
        0
    """
    NEG = -1
    ZERO = 0
    POS = 1

    def __str__(self):
        return {self.NEG: "-", self.ZERO: "0", self.POS: "+"}[self]

    def __repr__(self):
        return str(self)

    def __neg__(self):
        return Sign(-self.value) if self.value != 0 else Sign.ZERO

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            if value == "+":
                return cls.POS
            if value == "-":
                return cls.NEG
            return cls.ZERO
        v = sign_symbolic(value)
        if v > 0:
            return cls.POS
        if v < 0:
            return cls.NEG
        return cls.ZERO


class _Chirotope:
    def __init__(self, faces: set[SignVector], rank: int, ground_set_size: int):
        self.rank = rank
        self.ground_set_size = ground_set_size
        self._chirotope_dict: dict[tuple[int], Sign] = {}
        self._faces_dict: dict[tuple[int], SignVector] = {}
        self._set_faces_dict(faces)

    def _set_faces_dict(self, faces: set[SignVector]) -> None:
        r"""
        Set up a dictionary mapping index sets to faces.

        .. NOTE::

        TESTS::

            sage: from sign_vectors import *
            sage: from sign_vectors.chirotope import *
            sage: c = ChirotopeFromCircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
            sage: c._faces_dict
            {(0, 1, 2): (00+0), (0, 1, 3): (000+), (0, 2, 3): (000+), (1, 2, 3): (000+)}

        Note that the corresponding chirotopes are ``+00000``.
        Hence, ``(0, 2, 3)`` and ``(1, 2, 3)`` produce the zero face contrary to what is stored.
        However, this doesn't matter for the algorithm since these index sets connect two zero minors.

        We do the same for the cocircuits::

            sage: c = ChirotopeFromCocircuits([sign_vector("00+0"), sign_vector("000+")], 2, 4)
            sage: c._faces_dict
            {(0,): (000+), (1,): (000+), (2,): (000+), (3,): (00+0)}
        """
        checked_supports = set()
        for face in faces:
            support = frozenset(face.support())
            if support in checked_supports:
                continue
            checked_supports.add(support)
            for indices in self._corresponding_indices_for_face(support):
                self._faces_dict[indices] = face

    def _set_entry(self, rset: tuple[int], value: Sign) -> None:
        self._chirotope_dict[rset] = value

    def _has_entry(self, rset: list[int]) -> bool:
        return tuple(rset) in self._chirotope_dict

    def entry(self, rset: list[int]) -> Sign:
        r"""Return the chirotope entry given by ``indices``."""
        if self._chirotope_dict.get(tuple(rset)) is None:
            self._set_entries()
        return self._chirotope_dict.get(tuple(rset))

    def entries(self) -> list[Sign]:
        r"""Return all chirotope entries in lexicographic order."""
        return [self.entry(rset) for rset in Combinations(self.ground_set_size, self.rank)]

    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        r"""
        Generate all possible index sets that correspond to the face given by its support.
        
        For instance, ``{0, 1, 2}`` corresponds to the circuit ``(+++0)``.
        Also ``{0, 1, 2}`` would also correspond to the circuit ``(++00)``.

        On the other hand, ``{3}`` would correspond to the cocircuit ``(+++0)``.
        """
        raise NotImplementedError

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        raise NotImplementedError

    def _set_entries(self) -> None:
        self._set_other_nonzero_entries_from(self._set_zero_entries_and_first_nonzero())

    def _set_zero_entries_and_first_nonzero(self) -> tuple[int]:
        r"""
        Set all zero entries of the chirotope using (co)circuits and set the first nonzero entry to ``+``.
        """
        nonzero_rset = None
        for rset in Combinations(self.ground_set_size, self.rank):
            if self._is_entry_zero(rset):
                self._set_entry(tuple(rset), Sign.ZERO)
            elif nonzero_rset is None:
                nonzero_rset = tuple(rset)
                self._set_entry(nonzero_rset, Sign.POS)
        return nonzero_rset

    def _is_entry_zero(self, rset: list[int]) -> bool:
        r"""Return whether the entry of the chirotope given by ``rset`` is zero using (co)circuits."""
        raise NotImplementedError

    def _get_adjacent_rsets(self, rset: tuple[int]) -> Iterator[tuple[int]]:
        if self.rank == 0:
            return
        complement = [i for i in range(self.ground_set_size) if i not in rset]
        for subset in Combinations(rset, self.rank - 1):
            for i in complement:
                yield tuple(sorted(subset + [i]))

    def _set_other_nonzero_entries_from(self, rset: tuple[int]) -> None:
        r"""
        Recursively set all missing nonzero entries of the chirotope using the Grassmann-Plücker relations.
        """
        value = self.entry(rset)
        for adjacent_rset in self._get_adjacent_rsets(rset):
            if self._has_entry(adjacent_rset):
                continue

            face_indices = self._connecting_face_indices(rset, adjacent_rset)
            face = self._faces_dict[face_indices]

            # sort because when iterating over {7, 8}, we get 8 first
            i, j = sorted(set(rset).symmetric_difference(adjacent_rset))
            new_value = Sign(value * face[i] * face[j])
            if sum(i <= k < j for k in face_indices) & 1:
                new_value = -new_value
            self._set_entry(adjacent_rset, new_value)
            self._set_other_nonzero_entries_from(adjacent_rset)


class ChirotopeFromCircuits(_Chirotope):
    r"""A chirotope constructed from its circuits."""
    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        complement = set(range(self.ground_set_size)) - support
        for indices in Combinations(complement, self.rank + 1 - len(support)):
            yield tuple(sorted(support.union(indices)))

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        return tuple(sorted(set(rset1).union(rset2)))

    def _is_entry_zero(self, rset: list[int]) -> bool:
        for i in range(self.ground_set_size):
            if i in rset:
                continue
            face = self._faces_dict[tuple(sorted(rset + [i]))]
            if set(face.support()).issubset(rset):
                return True
        return False


class ChirotopeFromCocircuits(_Chirotope):
    r"""A chirotope constructed from its cocircuits."""
    def _corresponding_indices_for_face(self, support: set[int]) -> Iterator[tuple[int]]:
        complement = set(range(self.ground_set_size)) - support
        for indices in Combinations(complement, self.rank - 1):
            yield tuple(sorted(indices))

    def _connecting_face_indices(self, rset1: tuple[int], rset2: tuple[int]) -> tuple[int]:
        return tuple(sorted(set(rset1).intersection(rset2)))

    def _is_entry_zero(self, rset: list[int]) -> bool:
        for i in rset:
            face = self._faces_dict[tuple(sorted(set(rset) - {i}))]
            if set(rset).issubset(face.zero_support()):
                return True
        return False

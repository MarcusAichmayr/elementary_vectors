r"""Oriented matroids"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from typing import Iterator
from random import choice

from sage.combinat.combination import Combinations
from sage.combinat.posets.lattices import LatticePoset
from sage.structure.sage_object import SageObject

from sign_vectors import SignVector, sign_vector, zero_sign_vector
from .chirotope import Sign, ChirotopeFromCircuits, ChirotopeFromCocircuits
from .utility import classes_same_support, parallel_classes


class OrientedMatroid(SageObject):
    r"""Class representing an oriented matroid.

    This class provides methods to work with oriented matroids, including
    computing their faces, cocircuits, and topes.

    INPUT:

    - ``matrix``: A matrix representing the oriented matroid.
    - ``rank``: The rank of the oriented matroid. Only required if ``matrix`` is not provided.
    - ``ground_set_size``: The length of the elements in the oriented matroid. Only required if ``matrix`` is not provided.

    EXAMPLES:

    We define an oriented matroid from a matrix::

        sage: from sign_vectors import *
        sage: M = matrix([[1, 2, 0, 0, 0], [0, 1, 2, 0, 0], [0, 0, 0, 1, 1]])
        sage: om = OrientedMatroid(M)
        sage: om
        Oriented matroid of dimension 2 with elements of size 5.
        sage: om.ground_set
        {0, 1, 2, 3, 4}
        sage: om.dimension
        2
        sage: om.rank
        3

    We compute the cocircuits and topes::

        sage: om.cocircuits()
        {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)}
        sage: om.topes()
        {(-++++),
         (---++),
         (+--++),
         (++-++),
         (+----),
         (--+++),
         (+++++),
         (-----),
         (++---),
         (--+--),
         (-++--),
         (+++--)}
        sage: om.elements()
        {(00000),
         (+0---),
         (0----),
         (++-00),
         (--+00),
         (0++--),
         (000++),
         (++-++),
         (--+++),
         (0++++),
         (++---),
         (+--00),
         (--+--),
         (+0-++),
         (+----),
         (++000),
         (++0++),
         (++0--),
         (---00),
         (---++),
         (-0+00),
         (0++00),
         (-++00),
         (+--++),
         (+++00),
         (+++++),
         (000--),
         (-----),
         (-0+--),
         (-++--),
         (+++--),
         (-0+++),
         (-++++),
         (--000),
         (+0-00),
         (0--00),
         (0--++),
         (--0++),
         (--0--)}

    Next, we compute the faces of given dimensions::

        sage: om.faces(-1)
        {(00000)}
        sage: om.faces(0)
        {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)}
        sage: om.faces(1)
        {(+0---),
         (++0--),
         (---00),
         (0----),
         (-0+++),
         (++-00),
         (--+00),
         (0++--),
         (-++00),
         (+++00),
         (0--++),
         (0++++),
         (-0+--),
         (--0++),
         (++0++),
         (+--00),
         (--0--),
         (+0-++)}
        sage: om.faces(2)
        {(-++++),
         (---++),
         (+--++),
         (++-++),
         (+----),
         (--+++),
         (+++++),
         (-----),
         (++---),
         (--+--),
         (-++--),
         (+++--)}
        sage: om.faces()
        [{(00000)},
         {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)},
         {(+0---),
          (++0--),
          (---00),
          (0----),
          (-0+++),
          (++-00),
          (--+00),
          (0++--),
          (-++00),
          (+++00),
          (0--++),
          (0++++),
          (-0+--),
          (--0++),
          (++0++),
          (+--00),
          (--0--),
          (+0-++)},
         {(-++++),
          (---++),
          (+--++),
          (++-++),
          (+----),
          (--+++),
          (+++++),
          (-----),
          (++---),
          (--+--),
          (-++--),
          (+++--)}]

    We count the number of faces::

        sage: om.f_vector()
        [1, 8, 18, 12]
        sage: om.num_faces()
        39

    The dual oriented matroid corresponds to the kernel matrix.
    It is represented by the circuits and vectors::

        sage: om.circuits()
        {(000-+), (-+-00), (000+-), (+-+00)}
        sage: om.vectors()
        {(00000),
         (000+-),
         (-+-+-),
         (000-+),
         (+-+00),
         (+-+-+),
         (-+-00),
         (-+--+),
         (+-++-)}

    We compute the chirotope::

        sage: om.chirotope()
        [0, +, +, +, +, 0, +, +, 0, 0]
        sage: om.chirotope_string()
        '0++++0++00'
        sage: om.chirotope_entry([0, 1, 2])
        0
        sage: om.cocircuit([0, 1])
        (000++)
        sage: om.circuit([0, 1, 2, 3])
        (+-+00)

    Now, we compute the dual oriented matroid::

        sage: om_dual = om.dual()
        sage: om_dual
        Oriented matroid of dimension 1 with elements of size 5.
        sage: om_dual.cocircuits()
        {(000-+), (-+-00), (000+-), (+-+00)}
        sage: om.circuits()
        {(000-+), (-+-00), (000+-), (+-+00)}
        sage: om_dual.circuits()
        {(--000), (000++), (+0-00), (++000), (0--00), (0++00), (000--), (-0+00)}
        sage: om.cocircuits()
        {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)}
        sage: om_dual.topes()
        {(-+-+-), (+-+-+), (-+--+), (+-++-)}

    We compute the face lattice::

        sage: om.face_lattice()
        Finite lattice containing 40 elements

    TESTS:

    Trivial oriented matroid::

        sage: M = matrix(ZZ, 0, 3)
        sage: om = OrientedMatroid(M)
        sage: om.cocircuits()
        set()
        sage: om.circuits()
        {(-00), (00+), (00-), (0+0), (0-0), (+00)}
        sage: om.topes()
        {(000)}
        sage: om.faces()
        [{(000)}]
    """
    def __init__(self, matrix=None, rank: int = None, ground_set_size: int = None) -> None:
        if matrix is None:
            self._matrix = None
            self._rank = rank
            self._ground_set_size = ground_set_size
        else:
            try:
                self._matrix = matrix.matrix_from_rows(matrix.pivot_rows())
            except NotImplementedError as exc:
                if all(minor == 0 for minor in matrix.minors(matrix.nrows())):
                    raise ValueError("Provide a matrix with maximal rank.") from exc
                self._matrix = matrix
            self._rank, self._ground_set_size = self._matrix.dimensions()

        self._dimension = None
        self._chirotope_dict = {}
        self._faces_by_dimension = {}
        self._loops = None

        self._connect_faces = True
        self._above = {}
        self._connected_with_lower_dimension = set()  # faces of this dimensions are already connected with faces below

    def _repr_(self) -> str:
        return f"Oriented matroid of dimension {self.dimension} with elements of size {self.ground_set_size}."

    @classmethod
    def from_chirotope(cls, chirotope: list[int] | str, rank: int, ground_set_size: int) -> "OrientedMatroid":
        r"""
        Create an oriented matroid from a chirotope.

        INPUT:

        - ``chirotope`` -- a list of integers or a string of ``+``, ``-``, ``0``.

        OUTPUT:

        - An instance of :class:`OrientedMatroid`.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: om = OrientedMatroid.from_chirotope([1, 2, 3, 4, 6, 0], 2, 4)
            sage: om
            Oriented matroid of dimension 1 with elements of size 4.
            sage: om.faces()
            [{(0000)},
             {(+0--), (-0++), (--00), (++00), (0+++), (0---)},
             {(-+++), (+---), (--++), (++++), (----), (++--)}]

        ::

            sage: om = OrientedMatroid.from_chirotope([0, 0, -1, -1, 0, 1, 1, 1, 1, 1], 2, 5)
            sage: om
            Oriented matroid of dimension 1 with elements of size 5.
            sage: om.faces()
            [{(00000)},
             {(-+++0), (+---0), (000++), (-++0-), (+--0+), (000--)},
             {(-+++-), (-++++), (+--++), (+----), (-++--), (+---+)}]

        The chirotope can also be specified as a string::

            sage: OrientedMatroid.from_chirotope("00--0+++++", 2, 5)
            Oriented matroid of dimension 1 with elements of size 5.
        """
        om = cls(rank=rank, ground_set_size=ground_set_size)
        for (indices, value) in zip(Combinations(ground_set_size, rank), chirotope):
            om._set_chirotope_entry(tuple(indices), Sign(value))
        return om

    @classmethod
    def from_cocircuits(cls, cocircuits: set[SignVector | str], rank: int = None, ground_set_size: int = None) -> "OrientedMatroid":
        r"""
        Create an oriented matroid from cocircuits.

        INPUT:

        - ``cocircuits`` -- a set (or list) of sign vectors or strings representing the cocircuits.

        - ``rank`` -- the rank of the oriented matroid (optional).

        - ``ground_set_size`` -- the size of the ground set (optional).

        EXAMPLES::

            sage: from sign_vectors import *
            sage: om = OrientedMatroid.from_cocircuits({"0+", "+0"})
            sage: om
            Oriented matroid of dimension 1 with elements of size 2.
            sage: om.faces()
            [{(00)}, {(+0), (0-), (0+), (-0)}, {(--), (+-), (-+), (++)}]

        Specify `rank` to speed up computations::

            sage: om = OrientedMatroid.from_cocircuits({"0+", "+0"}, rank=2)
            sage: om
            Oriented matroid of dimension 1 with elements of size 2.
            sage: om.faces()
            [{(00)}, {(+0), (0-), (0+), (-0)}, {(--), (+-), (-+), (++)}]

        You can also use sign vector objects as input::

            sage: om = OrientedMatroid.from_cocircuits({sign_vector("0+"), sign_vector("+0")})
            sage: om
            Oriented matroid of dimension 1 with elements of size 2.
            sage: om.faces()
            [{(00)}, {(+0), (0-), (0+), (-0)}, {(--), (+-), (-+), (++)}]

        TESTS::

            sage: om = OrientedMatroid.from_cocircuits(["0++", "+00"])
            sage: om
            Oriented matroid of dimension 1 with elements of size 3.
            sage: sorted(om.faces(0), key=lambda X: hash(X))
            [(0++), (+00), (-00), (0--)]
            sage: om.faces(1)
            {(---), (+--), (+++), (-++)}

        ::

            sage: om = OrientedMatroid.from_cocircuits([], ground_set_size=3)
            sage: om
            Oriented matroid of dimension -1 with elements of size 3.
            sage: om.faces()
            [{(000)}]
            sage: om.circuits()
            {(-00), (00+), (00-), (0+0), (0-0), (+00)}
            sage: om.chirotope()
            [+]
            sage: om.dual()
            Oriented matroid of dimension 2 with elements of size 3.
            sage: om.dual().cocircuits()
            {(-00), (00+), (00-), (0+0), (0-0), (+00)}
        """
        def create_cocircuits(iterable):
            for element in iterable:
                yield sign_vector(element)
                yield -sign_vector(element)

        om = cls()
        if ground_set_size is None:
            try:
                om._ground_set_size = len(next(iter(cocircuits)))
            except StopIteration as e:
                raise ValueError("Could not determine 'ground_set_size'.") from e
        else:
            om._ground_set_size = ground_set_size
        if rank is None:
            if len(cocircuits) == 0:
                om._rank = 0
        else:
            om._rank = rank
        om._set_cocircuits(set(create_cocircuits(cocircuits)))
        return om

    @classmethod
    def from_circuits(cls, circuits: set[SignVector | str], rank: int, ground_set_size: int = None) -> "OrientedMatroid":
        r"""
        Create an oriented matroid from circuits.

        INPUT:

        - ``circuits`` -- a set (or list) of sign vectors or strings representing the circuits.

        - ``rank`` -- the rank of the oriented matroid.

        - ``ground_set_size`` -- the size of the ground set.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: om = OrientedMatroid.from_circuits({"0+0", "+00"}, 1)
            sage: om
            Oriented matroid of dimension 0 with elements of size 3.
            sage: om.chirotope()
            [0, 0, +]
            sage: om.faces()
            [{(000)}, {(00+), (00-)}]
        """
        om = cls()
        if ground_set_size is None:
            try:
                om._ground_set_size = len(next(iter(circuits)))
            except StopIteration as e:
                raise ValueError("Could not determine 'ground_set_size'.") from e
        else:
            om._ground_set_size = ground_set_size
        om._rank = rank
        om._set_chirotope_entries_from_circuits([sign_vector(element) for element in circuits])
        return om

    @classmethod
    def from_topes(cls, topes: list[SignVector | str] | set[SignVector | str]) -> "OrientedMatroid":
        r"""
        Create an oriented matroid from a tope set.

        INPUT:

        - ``topes`` -- a list or set of sign vectors or strings representing the topes.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: om = OrientedMatroid.from_topes({"+++"})
            sage: om
            Oriented matroid of dimension 0 with elements of size 3.
            sage: om.faces()
            [{(000)}, {(---), (+++)}]

        ::

            sage: om = OrientedMatroid.from_topes([[-1, 1, -1, 1], [-1, -1, -1, 1], [-1, 1, -1, -1], [-1, -1, 1, 1]])
            sage: om
            Oriented matroid of dimension 1 with elements of size 4.
            sage: om.faces()
            [{(0000)},
             {(++0-), (+0+-), (--0+), (0-++), (0+--), (-+-0), (-0-+), (+-+0)},
             {(-+-+), (---+), (-+--), (--++), (+-+-), (++--), (+-++), (+++-)}]
        """
        def create_topes(iterable):
            for element in iterable:
                yield sign_vector(element)
                yield -sign_vector(element)

        topes = set(create_topes(topes))

        om = cls()
        om._ground_set_size = next(iter(topes)).length()
        om._set_faces_from_topes(topes)
        om._rank = max(om._faces_by_dimension) + 1
        return om

    @property
    def matrix(self):
        r"""The underlying matrix if available."""
        if self._matrix is None:
            raise ValueError("There is no matrix associated with this oriented matroid.")
        return self._matrix

    @property
    def rank(self) -> int:
        r"""The rank of this oriented matroid."""
        if self._rank is None:
            self._rank = self._rank_from_cocircuits(self.cocircuits())
        return self._rank

    def _rank_from_cocircuits(self, cocircuits: set[SignVector]) -> int:
        self._set_faces_from_topes(self._topes_from_cocircuits(cocircuits))
        return max(self._faces_by_dimension) + 1

    @property
    def dimension(self) -> int:
        r"""The dimension of this oriented matroid."""
        if self._dimension is None:
            self._dimension = self.rank - 1
        return self._dimension

    @property
    def ground_set_size(self) -> int:
        r"""The size of the ground set of this oriented matroid."""
        if self._ground_set_size is None:
            if not self._faces_by_dimension:
                raise ValueError("Couldn't determine ground set size.")
            for faces in self._faces_by_dimension.values():
                if faces:
                    self._ground_set_size = next(iter(faces)).length()
                    break
        return self._ground_set_size

    @property
    def ground_set(self) -> set[int]:
        r"""The ground set of this oriented matroid."""
        return set(range(self.ground_set_size))

    def loops(self) -> list[int]:
        r"""
        Compute the loops of this oriented matroid.

        A loop is a component where every face is zero.

        .. NOTE::

            The result is cached.
        """
        if self._loops is None:
            self._loops = [
                e for e in range(self.ground_set_size)
                if all(element[e] == 0 for element in self.cocircuits())
            ]
        return self._loops

    def parallel_classes(self) -> list[set[int]]:
        r"""Compute the parallel classes of this oriented matroid."""
        return parallel_classes(self.cocircuits(), self.ground_set_size)

    def chirotope_entry(self, indices: list[int]) -> Sign:
        r"""
        Compute the chirotope for the given indices.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The chirotope value as a ``Sign``.

        .. SEEALSO::

            - :meth:`chirotopes`
            - :class:`Sign`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotope_entry([1, 2])
            +
        """
        indices = tuple(indices)
        if self._chirotope_dict.get(indices) is None:
            try:
                if self._matrix is None:
                    self._set_chirotope_entries_from_cocircuits()
                else:
                    self._set_chirotope_entry(
                        indices,
                        Sign(self.matrix.matrix_from_columns(indices).det())
                    )
            except ValueError as e:
                if len(indices) != self.rank:
                    raise ValueError(f"Indices {indices} should have size {self.rank} and not {len(indices)}.") from e
                raise e
        return self._chirotope_dict.get(indices)

    def _set_chirotope_entry(self, indices: tuple[int], value: Sign) -> None:
        r"""
        Set the chirotope for the given indices.
        
        INPUT:

        - ``indices`` -- a sorted tuple of integers

        - ``value`` -- the chirotope value as a ``Sign``.
        """
        self._chirotope_dict[indices] = value

    def _set_chirotope_entries_from_cocircuits(self) -> None:
        chirotope = ChirotopeFromCocircuits(self.cocircuits(), self.rank, self.ground_set_size)
        for rset in Combinations(self.ground_set_size, self.rank):
            value = chirotope.entry(rset)
            self._set_chirotope_entry(tuple(rset), value)

    def _set_chirotope_entries_from_circuits(self, circuits: set[SignVector]) -> None:
        chirotope = ChirotopeFromCircuits(circuits, self.rank, self.ground_set_size)
        for rset in Combinations(self.ground_set_size, self.rank):
            value = chirotope.entry(rset)
            self._set_chirotope_entry(tuple(rset), value)

    def _compute_chirotope_entries(self) -> None:
        for indices in Combinations(self.ground_set_size, self.rank):
            self.chirotope_entry(indices)

    def chirotope(self) -> list[Sign]:
        r"""
        Compute all chirotope entries of the oriented matroid.

        OUTPUT:

        - A list of chirotope entries as ``Sign`` values.

        .. SEEALSO::

            - :meth:`chirotope_entry`
            - :class:`Sign`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotope()
            [+, +, +, +, +, 0]
        """
        return [self.chirotope_entry(indices) for indices in Combinations(self.ground_set_size, self.rank)]

    def chirotope_string(self) -> str:
        r"""Represent the chirotope as a string."""
        return "".join(str(value) for value in self.chirotope())

    def dual(self) -> "OrientedMatroid":
        r"""
        Return the dual oriented matroid.

        .. NOTE::

            The dual is determined from the chirotope.
        """
        self._compute_chirotope_entries()
        om = OrientedMatroid(rank=self.ground_set_size - self.rank, ground_set_size=self.ground_set_size)
        for indices, value in self._chirotope_dict.items():
            indices_set = set(indices)
            complement = tuple(i for i in range(self.ground_set_size) if i not in indices_set)
            inversions = sum(i < j for i in indices for j in complement)
            om._set_chirotope_entry(complement, Sign(-value if inversions & 1 else value)) # check last bit
        return om

    def _set_zero_face(self) -> None:
        r"""Set the zero face of the oriented matroid."""
        self._set_faces(-1, {zero_sign_vector(self.ground_set_size)})

    def cocircuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute a cocircuit for the given indices.

        Cocircuits correspond to the elements in the row space of the matrix.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The cocircuit as a sign vector of this oriented matroid.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`circuit`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.cocircuit([0])
            (0---)
        """
        element = [0] * self.ground_set_size
        pos = 0
        for i in range(self.ground_set_size):
            if i in indices:
                pos += 1
                continue
            value = self.chirotope_entry(sorted(indices + [i]))
            if value != 0:
                # check whether pos is even or odd
                element[i] = -value if (pos & 1) else value
        result = sign_vector(element)
        if result == 0:
            raise ValueError("Computed zero cocircuit.")
        return result

    def circuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute a circuit for the given indices.

        Circuits correspond to the elements in the kernel of the matrix.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The circuit as a sign vector of this oriented matroid.

        .. SEEALSO::

            - :meth:`circuits`
            - :meth:`cocircuit`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.circuit([0, 1, 2])
            (+-+0)
        """
        element = [0] * self.ground_set_size
        for pos in range(self.rank + 1):
            rset = indices.copy()
            i = rset.pop(pos)
            value = self.chirotope_entry(rset)
            if value != 0:
                # check whether pos is even or odd
                element[i] = -value if (pos & 1) else value
        result = sign_vector(element)
        if result == 0:
            raise ValueError("Computed zero circuit.")
        return result

    def _cocircuit_generator(self) -> Iterator[SignVector]:
        if self.rank == 0:
            return
        for indices in Combinations(self.ground_set_size, self.rank - 1):
            try:
                yield self.cocircuit(indices)
                yield -self.cocircuit(indices)
            except ValueError:
                continue

    def _circuit_generator(self) -> Iterator[SignVector]:
        for indices in Combinations(self.ground_set_size, self.rank + 1):
            try:
                yield self.circuit(indices)
                yield -self.circuit(indices)
            except ValueError:
                continue

    def cocircuits(self) -> set[SignVector]:
        r"""
        Compute the cocircuits of the oriented matroid.

        The cocircuits are the support-minimal sign vectors corresponding to the row space.

        OUTPUT:

        - A set of cocircuits as sign vectors of this oriented matroid.

        .. NOTE::

            The result is cached.

        .. SEEALSO::

            - :meth:`cocircuit`
            - :meth:`circuits`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.cocircuits()
            {(+0--), (-0++), (--00), (++00), (0+++), (0---)}
        """
        if self._rank is not None:
            if self.dimension == -1:
                return set()
            if 0 not in self._faces_by_dimension:
                self._set_cocircuits(set(self._cocircuit_generator()))
        return self._faces_by_dimension[0]

    def _set_cocircuits(self, cocircuits: set[SignVector]) -> None:
        if len(cocircuits) == 0:
            return
        self._set_faces(0, cocircuits)

    def circuits(self) -> set[SignVector]:
        r"""
        Compute the circuits of the oriented matroid.

        OUTPUT:

        - A set of circuits as sign vectors of this oriented matroid.

        .. NOTE::

            The circuits are the support-minimal sign vectors corresponding to the kernel.

        .. SEEALSO::

            - :meth:`circuit`
            - :meth:`cocircuits`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.circuits()
            {(00-+), (+-0+), (-+0-), (00+-), (-+-0), (+-+0)}
        """
        return set(self._circuit_generator())

    def an_element(self) -> SignVector:
        r"""
        Compute a random element of the oriented matroid.

        .. SEEALSO::

            - :meth:`elements`
        """
        return choice(list(self.elements()))

    def covectors(self) -> set[SignVector]:
        r"""
        Compute the covectors of the oriented matroid.

        .. NOTE::

            The covectors are all sign vectors corresponding to the row space.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`topes`
            - :meth:`faces`
        """
        if self._topes_computed():
            return set().union(*self._all_faces())
        return self._faces_from_vertices(self.cocircuits())

    def elements(self) -> set[SignVector]:
        r"""
        Compute the covectors of the oriented matroid.

        .. SEEALSO::

            - :meth:`covectors`
        """
        return self.covectors()

    def vectors(self) -> set[SignVector]:
        r"""
        Compute the vectors of the oriented matroid.

        .. NOTE::

            The vectors are all sign vectors corresponding to the kernel.

        .. SEEALSO::

            - :meth:`circuits`
        """
        return self._faces_from_vertices(self.circuits())

    def topes(self) -> set[SignVector]:
        r"""
        Compute the topes of the oriented matroid.

        The topes are the support-maximal sign vectors corresponding to the row space.

        .. SEEALSO::

            - :meth:`cocircuits`
        """
        if not self._topes_computed():
            if self.dimension == -1:
                self._set_zero_face()
            else:
                self._set_faces(self.dimension, self._topes_from_cocircuits(self.cocircuits()))
        return self._faces_by_dimension[self.dimension]

    def faces(self, dimension: int = None) -> set[SignVector] | list[set[SignVector]]:
        r"""
        Compute the faces of the same level of the oriented matroid.

        OUTPUT:
        If ``dimension`` is not specified, return all faces separated by dimension.

        .. NOTE::

            The result is cached.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`topes`
        """
        if dimension is None:
            return self._all_faces()
        if dimension < -1 or dimension > self.dimension:
            raise ValueError(f"Dimension should be between -1 and {self.dimension}. Got {dimension}.")
        if dimension in self._faces_by_dimension:
            return self._faces_by_dimension[dimension]
        if dimension == 0:
            return self.cocircuits()

        self.topes()
        current_dimension = self.dimension
        while current_dimension > dimension:
            self._set_lower_faces(current_dimension)
            current_dimension -= 1
        return self._faces_by_dimension[dimension]

    def _set_faces(self, dimension: int, faces: set[SignVector]) -> None:
        if faces == set():
            raise ValueError("Attempting to set empty face set.")
        self._faces_by_dimension[dimension] = faces

    def _all_faces(self) -> list[set[SignVector]]:
        r"""
        Return all faces of the oriented matroid separated by dimension.

        .. SEEALSO::

            - :meth:`faces`
        """
        return [self.faces(d) for d in range(-1, self.dimension + 1)]

    def num_faces(self) -> int:
        r"""
        Return the total number of faces (covectors) of the oriented matroid.

        .. SEEALSO::

            - :meth:`f_vector`
            - :meth:`faces`
        """
        return sum(len(faces) for faces in self._all_faces())

    def f_vector(self) -> list[int]:
        r"""
        Compute the f-vector of the oriented matroid.

        The f-vector (face vector) is a list where the ``i``-th entry is the number of faces of dimension ``i``.

        .. SEEALSO::

            - :meth:`num_faces`
            - :meth:`faces`
        """
        return [len(faces) for faces in self._all_faces()]

    def _topes_from_cocircuits(self, cocircuits: set[SignVector]) -> set[SignVector]:
        r"""
        Compute the topes from the cocircuits.

        OUTPUT:
        A set of topes of the oriented matroid.

        ALGORITHM:
        This function is based on an algorithm in [Fin01]_.
        """
        covectors = {zero_sign_vector(self.ground_set_size)}
        covectors_new = {zero_sign_vector(self.ground_set_size)}
        topes = set()

        while covectors_new:
            element1 = covectors_new.pop()
            for element2 in cocircuits:
                if element2 <= element1:
                    continue
                new_element = element2.compose(element1)
                if new_element not in covectors:
                    covectors.add(new_element)
                    if new_element.zero_support() == self.loops():
                        topes.add(new_element)
                    else:
                        covectors_new.add(new_element)
        return topes

    def _topes_computed(self) -> bool:
        return self.dimension in self._faces_by_dimension

    def _faces_from_vertices(self, vertices: set[SignVector]) -> set[SignVector]:
        r"""
        Compute the covectors from the cocircuits (vertices).

        OUTPUT:

        - a set of all covectors of the oriented matroid.

        ALGORITHM:

        This function is based on ``CovectorsFromCocircuits`` from [Fin01]_.

        .. [Fin01] Finschi, L.:
        „A graph theoretical approach for reconstruction and generation of oriented matroids“.
        PhD thesis. Zurich: ETH Zurich, 2001. doi: 10.3929/ethz-a-004255224.
        """
        covectors = {zero_sign_vector(self.ground_set_size)}
        covectors_new = {zero_sign_vector(self.ground_set_size)}
        while covectors_new:
            element1 = covectors_new.pop()
            for element2 in vertices:
                if element2 <= element1:
                    continue
                new_element = element2.compose(element1)
                if new_element not in covectors:
                    covectors.add(new_element)
                    covectors_new.add(new_element)
        return covectors

    def _set_faces_from_topes(self, topes: set[SignVector]) -> None:
        r"""
        Use the topes to set all faces of the oriented matroid.

        INPUT:

        - ``topes`` -- a set of topes of the oriented matroid

        ALGORITHM:

        This function is based on an algorithm in [FST91]_.
        See also [Fin01]_.
        """
        all_faces = [set(topes)]

        current = 0
        while len(all_faces[current]) > 1:
            all_faces.append(self._lower_faces(all_faces[current], self._connect_faces))
            current += 1

        dimension = -1
        while all_faces:
            self._set_faces(dimension, all_faces.pop())
            dimension += 1

    def _lower_faces(self, faces: set[SignVector], connect_faces: bool) -> set[SignVector]:
        r"""
        Compute the faces of lower dimension.

        INPUT:

        - ``dimension`` -- the dimension ``i``

        OUTPUT:
        Return a set of faces of dimension ``i - 1`` of the oriented matroid.

        ALGORITHM:

        This function is based on an algorithm in [FST91]_.
        See also [Fin01]_.

        .. [FST91] Fukuda, K., Saito, S., and Tamura, A.:
        "Combinatorial face enumeration in arrangements and oriented matroids".
        In: Discrete Applied Mathematics 31.2 (1991), pp. 141-149.
        doi: 10.1016/0166-218X(91)90066-6.

        .. SEEALSO::

            - :meth:`faces`
        """
        output = set()
        for same_support_faces in classes_same_support(faces):
            p_classes = parallel_classes(same_support_faces, self.ground_set_size)
            while same_support_faces:
                face = same_support_faces.pop()
                for parallel_class in p_classes:
                    if all(face[i] == 0 for i in parallel_class):
                        continue
                    flipped_face = face.flip_signs(parallel_class)
                    if flipped_face in same_support_faces:
                        lower_face = face.set_to_zero(parallel_class)
                        output.add(lower_face)
                        output.add(-lower_face)
                        if connect_faces:
                            self._connect(lower_face, face)
                            self._connect(-lower_face, -face)
                            self._connect(lower_face, flipped_face)
                            self._connect(-lower_face, -flipped_face)
                same_support_faces.remove(-face)
        return output

    def _set_lower_faces(self, dimension: int) -> None:
        r"""Set faces for one lower dimension."""
        if not self._faces_by_dimension.get(dimension):
            raise ValueError(f"Dimension {dimension} is not available. Available dimensions: {sorted(self._faces_by_dimension.keys())}.")
        if dimension - 1 in self._faces_by_dimension:
            return
        if dimension == 0:
            self._set_zero_face()
            return
        connect_faces = self._connect_faces and dimension not in self._connected_with_lower_dimension
        self._set_faces(dimension - 1, self._lower_faces(self._faces_by_dimension[dimension], connect_faces))
        if connect_faces:
            self._connected_with_lower_dimension.add(dimension)

    def set_face_connections(self, value: bool = True) -> None:
        """
        Keep track of connections between faces.

        If this is set to true, the object will keep track of connections between faces while computing.
        This improves performance for plotting but slightly slows down other computations.
        """
        self._connect_faces = value

    def _connect(self, lower_face: SignVector, upper_face: SignVector):
        r"""Connect two faces."""
        if lower_face not in self._above:
            self._above[lower_face] = set()
        self._above[lower_face].add(upper_face)

    def _connect_below(self, dimension: int) -> None:
        if dimension in self._connected_with_lower_dimension:
            return
        if dimension - 1 not in self._faces_by_dimension:
            raise ValueError(f"Trying to connect faces of dimension {dimension - 1} and {dimension}, but dimension {dimension - 1} is not available.")
        if dimension == self.dimension + 1:
            for tope in self.topes():
                self._connect(tope, 1)
            return
        if dimension not in self._faces_by_dimension:
            raise ValueError(f"Trying to connect faces of dimension {dimension - 1} and {dimension}, but dimension {dimension} is not available.")
        if dimension == 0:
            zero = zero_sign_vector(self.ground_set_size)
            for face in self.cocircuits():
                self._connect(zero, face)
        elif dimension == 1:
            for face in self._faces_by_dimension[1]:
                connection_count = 0
                for cocircuit in self.cocircuits():
                    if cocircuit <= face:
                        self._connect(cocircuit, face)
                        connection_count += 1
                        if connection_count == 2: # diamond property
                            break
        else:
            for face in self._faces_by_dimension[dimension]:
                for lower_face in self._faces_by_dimension[dimension - 1]:
                    if lower_face <= face:
                        self._connect(lower_face, face)
        self._connected_with_lower_dimension.add(dimension)

    def _connect_all(self):
        r"""
        Connect all faces of the oriented matroid.

        TESTS::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
            sage: om = OrientedMatroid(M)
            sage: om.faces(2)
            {(---0), (++-0), (--+0), (-++0), (+++0), (+--0), (-+-0), (+-+0)}
            sage: om._above
            {}
            sage: om.faces(1)
            {(0+-0),
             (+-00),
             (-+00),
             (--00),
             (0-+0),
             (0--0),
             (++00),
             (+0-0),
             (+0+0),
             (-0-0),
             (0++0),
             (-0+0)}
            sage: om._above
            {(0--0): {(---0), (+--0)},
             (0++0): {(-++0), (+++0)},
             (-0-0): {(---0), (-+-0)},
             (+0+0): {(+++0), (+-+0)},
             (--00): {(---0), (--+0)},
             (++00): {(+++0), (++-0)},
             (0+-0): {(-+-0), (++-0)},
             (0-+0): {(+-+0), (--+0)},
             (+0-0): {(+--0), (++-0)},
             (-0+0): {(-++0), (--+0)},
             (-+00): {(-++0), (-+-0)},
             (+-00): {(+--0), (+-+0)}}
            sage: om._connect_all()
            sage: om._above
            {(0000): {(-000), (00+0), (00-0), (0+00), (0-00), (+000)},
             (0--0): {(---0), (+--0)},
             (00+0): {(-0+0), (+0+0), (0-+0), (0++0)},
             (0++0): {(-++0), (+++0)},
             (-0-0): {(---0), (-+-0)},
             (+0+0): {(+++0), (+-+0)},
             (--00): {(---0), (--+0)},
             (0+00): {(0+-0), (0++0), (-+00), (++00)},
             (++00): {(+++0), (++-0)},
             (0+-0): {(-+-0), (++-0)},
             (0-+0): {(+-+0), (--+0)},
             (+0-0): {(+--0), (++-0)},
             (-0+0): {(-++0), (--+0)},
             (-+00): {(-++0), (-+-0)},
             (0-00): {(+-00), (0--0), (0-+0), (--00)},
             (+000): {(+-00), (+0+0), (+0-0), (++00)},
             (+-00): {(+--0), (+-+0)},
             (00-0): {(0+-0), (-0-0), (+0-0), (0--0)},
             (-000): {(--00), (-0+0), (-+00), (-0-0)},
             (---0): {1},
             (++-0): {1},
             (--+0): {1},
             (-++0): {1},
             (+++0): {1},
             (+--0): {1},
             (-+-0): {1},
             (+-+0): {1}}

        For type consistency, the top element ``1`` should not appear in the faces::

            sage: len(om._faces_by_dimension)
            4

        ::

            sage: om = OrientedMatroid(M)
            sage: om.set_face_connections(False)
            sage: om.faces()
            [{(0000)},
             {(-000), (00+0), (00-0), (0+00), (0-00), (+000)},
             {(0+-0),
              (+-00),
              (-+00),
              (--00),
              (0-+0),
              (0--0),
              (++00),
              (+0-0),
              (+0+0),
              (-0-0),
              (0++0),
              (-0+0)},
             {(---0), (++-0), (--+0), (-++0), (+++0), (+--0), (-+-0), (+-+0)}]
            sage: om._above
            {}
        """
        self._all_faces()
        for dimension in range(self.rank + 1):
            self._connect_below(dimension)

    def face_lattice(self) -> LatticePoset:
        r"""
        Return the face lattice of the oriented matroid.

        OUTPUT:

        - A lattice representing the face lattice.

        TESTS::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 1, 1]])
            sage: om = OrientedMatroid(M)
            sage: om.face_lattice()
            Finite lattice containing 4 elements

        ::

            sage: M = matrix(0, 3)
            sage: om = OrientedMatroid(M)
            sage: om.face_lattice()
            Finite lattice containing 2 elements

        ::

            sage: M = matrix([[1, 0], [0, 1]])
            sage: om = OrientedMatroid(M)
            sage: om.face_lattice()
            Finite lattice containing 10 elements
        """
        self.set_face_connections(True)
        self._connect_all()
        return LatticePoset(self._above)

    def plot(self, vertex_size: int = 600, figsize: int = None, aspect_ratio=None) -> None:
        r"""
        Plot the big face lattice of the oriented matroid.

        INPUT:

        - ``vertex_size`` -- the size of the vertices in the plot.

        - ``figsize`` -- the size of the figure.

        - ``aspect_ratio`` -- the aspect ratio of the plot.

        .. NOTE::

            Only works well for small length and dimension.
        """
        self.face_lattice().plot(
            vertex_size=vertex_size,
            element_color="white",
            vertex_shape="",
        ).show(figsize=figsize, aspect_ratio=aspect_ratio)

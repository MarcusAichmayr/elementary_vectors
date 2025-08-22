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

from collections.abc import Generator
from enum import IntEnum

from sage.combinat.combination import Combinations
from sage.combinat.posets.posets import Poset
from sage.structure.sage_object import SageObject

from sign_vectors import sign_symbolic, SignVector, sign_vector, zero_sign_vector

from .utility import classes_same_support, parallel_classes
from .functions import plot_sign_vectors


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


class OrientedMatroid(SageObject):
    r"""
    Class representing an oriented matroid.

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
        sage: om.covectors()
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
        sage: om.faces_all()
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

    Geometrically, the cocircuits are the vertices and the edges are the faces of dimension 1::

        sage: om.vertices()
        {(--000),
         (000++),
         (+0-00),
         (0--00),
         (++000),
         (0++00),
         (000--),
         (-0+00)}
        sage: om.edges()
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

    We compute the chirotopes::

        sage: om.chirotopes()
        [0, +, +, +, +, 0, +, +, 0, 0]
        sage: om.chirotope_string()
        '0++++0++00'
        sage: om.chirotope([0, 1, 2])
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
        sage: om.faces_all()
        [{(000)}]
    """
    def __init__(self, matrix=None, rank: int = None, element_length: int = None) -> None:
        if matrix is None:
            if rank is None or element_length is None:
                raise ValueError("Provide either a matrix or both rank and element_length.")
            self._matrix = None
            self.rank = rank
            self._element_length = element_length
        else:
            try:
                self._matrix = matrix.matrix_from_rows(matrix.pivot_rows())
            except NotImplementedError as exc:
                if all(minor == 0 for minor in matrix.minors(matrix.nrows())):
                    raise ValueError("Provide a matrix with maximal rank.") from exc
                self._matrix = matrix
            self.rank, self._element_length = self._matrix.dimensions()

        self.dimension = self.rank - 1

        self._chirotope_dict = {}
        self._faces_by_dimension = {-1: set([zero_sign_vector(self._element_length)])}
        self._loops = None

    def _repr_(self) -> str:
        return f"Oriented matroid of dimension {self.dimension} with elements of size {self._element_length}."

    @property
    def ground_set(self) -> set[int]:
        r"""The ground set of this oriented matroid."""
        return set(range(self._element_length))

    def chirotope(self, indices: list[int]) -> Sign:
        r"""
        Compute the chirotope for the given indices.

        Chirotopes are signs of determinants of maximal submatrices.

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
            sage: om.chirotope([1, 2])
            +
        """
        indices = tuple(indices)
        chirotope = self._chirotope_dict.get(indices)
        if chirotope is None:
            try:
                chirotope = Sign(self._matrix.matrix_from_columns(indices).det())
                self._chirotope_dict[indices] = chirotope
            except ValueError as e:
                raise ValueError(f"Indices {indices} should have size {self.rank} and not {len(indices)}.") from e
        return chirotope

    def set_chirotope(self, indices: list[int], value: Sign) -> None:
        r"""Set the chirotope for the given indices."""
        self._chirotope_dict[tuple(indices)] = Sign(value)

    def chirotopes(self) -> list[Sign]:
        r"""
        Compute all chirotopes of the oriented matroid.

        OUTPUT:

        - A list of chirotopes as ``Sign`` values.

        .. SEEALSO::

            - :meth:`chirotope`
            - :class:`Sign`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotopes()
            [+, +, +, +, +, 0]
        """
        return [self.chirotope(indices) for indices in Combinations(self._element_length, self.rank)]

    def chirotope_string(self) -> str:
        r"""Represent the chirotopes as a string."""
        return "".join(str(chirotope) for chirotope in self.chirotopes())

    def dual(self) -> "OrientedMatroid":
        r"""
        Return the dual oriented matroid.

        .. NOTE::

            The dual is determined from the chirotopes.
        """
        self.chirotopes() # compute all chirotopes
        om = OrientedMatroid(rank=self._element_length - self.rank, element_length=self._element_length)
        for indices, value in self._chirotope_dict.items():
            indices_set = set(indices)
            complement = tuple(i for i in range(self._element_length) if i not in indices_set)
            inversions = sum(i < j for i in indices for j in complement)
            om.set_chirotope(complement, -value if inversions & 1 else value) # check last bit
        return om

    def cocircuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute a cocircuit for the given indices.

        Cocircuits correspond to the elements in the row space of the matrix.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The cocircuit as a SignVector.

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
        element = [0] * self._element_length
        pos = 0
        for i in range(self._element_length):
            if i in indices:
                pos += 1
                continue
            chirotope = self.chirotope(sorted(indices + [i]))
            if chirotope != 0:
                # check oddness of last bit of pos
                element[i] = -chirotope if (pos & 1) else chirotope
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

        - The circuit as a SignVector.

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
        element = [0] * self._element_length
        for pos in range(self.rank + 1):
            indices_chirotope = indices.copy()
            i = indices_chirotope.pop(pos)
            chirotope = self.chirotope(indices_chirotope)
            if chirotope != 0:
                # check oddness of last bit of pos
                element[i] = -chirotope if (pos & 1) else chirotope
        result = sign_vector(element)
        if result == 0:
            raise ValueError("Computed zero circuit.")
        return result

    def _cocircuit_generator(self) -> Generator[SignVector]:
        r"""
        Compute the cocircuits of the oriented matroid.

        OUTPUT:

        - A generator of cocircuits as SignVectors.

        .. SEEALSO::

            - :meth:`cocircuits`
        """
        if self.rank == 0:
            return
        for indices in Combinations(self._element_length, self.rank - 1):
            try:
                yield self.cocircuit(indices)
                yield -self.cocircuit(indices)
            except ValueError:
                continue

    def cocircuits(self) -> set[SignVector]:
        r"""
        Compute the cocircuits of the oriented matroid.

        The cocircuits are the support-minimal sign vectors corresponding to the row space.

        OUTPUT:

        - A set of cocircuits as SignVectors.

        .. NOTE::

            The result is hashed.

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
        if self.dimension == -1:
            return set()
        if 0 not in self._faces_by_dimension:
            self._faces_by_dimension[0] = set(self._cocircuit_generator())
        return self._faces_by_dimension[0]

    def _circuit_generator(self) -> Generator[SignVector]:
        r"""
        Compute the circuits of the oriented matroid.

        OUTPUT:

        - A generator of circuits as SignVectors.

        .. SEEALSO::

            - :meth:`circuits`
        """
        for indices in Combinations(self._element_length, self.rank + 1):
            try:
                yield self.circuit(indices)
                yield -self.circuit(indices)
            except ValueError:
                continue

    def circuits(self) -> set[SignVector]:
        r"""
        Compute the circuits of the oriented matroid.

        OUTPUT:

        - A set of circuits as SignVectors.

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
            return set().union(*self.faces_all())
        return self.faces_from_vertices(self.cocircuits(), self._element_length)

    def vectors(self) -> set[SignVector]:
        r"""
        Compute the vectors of the oriented matroid.

        .. NOTE::

            The vectors are all sign vectors corresponding to the kernel.

        .. SEEALSO::

            - :meth:`circuits`
        """
        return self.faces_from_vertices(self.circuits(), self._element_length)

    def topes(self) -> set[SignVector]:
        r"""
        Compute the topes of the oriented matroid.

        The topes are the support-maximal sign vectors corresponding to the row space.

        .. SEEALSO::

            - :meth:`cocircuits`
        """
        if self._faces_by_dimension.get(self.dimension) is None:
            self._faces_by_dimension[self.dimension] = self._topes_from_cocircuits()
        return self._faces_by_dimension[self.dimension]

    def vertices(self) -> set[SignVector]:
        r"""
        Return the vertices (cocircuits) of the oriented matroid.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`edges`
            - :meth:`topes`
            - :meth:`faces`
        """
        return self.cocircuits()

    def edges(self) -> set[SignVector]:
        r"""
        Return the edges of the oriented matroid.

        Those are the elements of dimension 1.

        .. SEEALSO::

            - :meth:`vertices`
            - :meth:`topes`
            - :meth:`faces`
        """
        return self.faces(1)

    def faces(self, dimension: int) -> set[SignVector]:
        r"""
        Compute the faces of the same level of the oriented matroid.

        .. NOTE::

            The result is hashed.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`topes`
            - :meth:`faces_all`
        """
        if dimension < -1 or dimension > self.dimension:
            raise ValueError(f"Dimension should be between -1 and {self.dimension}. Got {dimension}.")
        if dimension in self._faces_by_dimension:
            return self._faces_by_dimension[dimension]
        if dimension == 0:
            return self.cocircuits()

        self.topes()
        current_dimension = self.dimension
        while current_dimension > dimension:
            self._compute_lower_faces(current_dimension)
            current_dimension -= 1
        return self._faces_by_dimension[dimension]

    def faces_all(self) -> list[set[SignVector]]:
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
            - :meth:`faces_all`
        """
        return sum(len(faces) for faces in self.faces_all())

    def f_vector(self) -> list[int]:
        r"""
        Compute the f-vector of the oriented matroid.

        The f-vector (face vector) is a list where the ``i``-th entry is the number of faces of dimension ``i``.

        .. SEEALSO::

            - :meth:`num_faces`
            - :meth:`faces`
        """
        return [len(faces) for faces in self.faces_all()]

    def _topes_computed(self) -> bool:
        return self.dimension in self._faces_by_dimension

    def _compute_lower_faces(self, dimension: int) -> None:
        if dimension - 1 not in self._faces_by_dimension:
            self._faces_by_dimension[dimension - 1] = self._lower_faces(dimension)

    def _lower_faces(self, dimension: int) -> set[SignVector]:
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
        if not self._faces_by_dimension.get(dimension):
            raise ValueError(f"Dimension {dimension} is not available. Available dimensions: {sorted(self._faces_by_dimension.keys())}.")
        if dimension == -1:
            return set()
        output = set()
        for same_support_faces in classes_same_support(self._faces_by_dimension[dimension]):
            p_classes = parallel_classes(same_support_faces, self._element_length)
            while same_support_faces:
                face = same_support_faces.pop()
                for parallel_class in p_classes:
                    if all(face[i] == 0 for i in parallel_class):
                        continue
                    if face.flip_signs(parallel_class) in same_support_faces:
                        lower_face = face.set_to_zero(parallel_class)
                        output.add(lower_face)
                        output.add(-lower_face)
                same_support_faces.remove(-face)
        return output

    def _topes_from_cocircuits(self) -> set[SignVector]:
        r"""
        Compute the topes from the cocircuits.

        OUTPUT:
        A set of topes of the oriented matroid.

        ALGORITHM:
        This function is based on an algorithm in [Fin01]_.
        """
        covectors = {zero_sign_vector(self._element_length)}
        covectors_new = {zero_sign_vector(self._element_length)}
        topes = set()

        while covectors_new:
            element1 = covectors_new.pop()
            for element2 in self.cocircuits():
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

    def loops(self) -> list[int]:
        r"""
        Compute the loops of this oriented matroid.

        A loop is a component where every face is zero.

        .. NOTE::

            The result is hashed.
        """
        if self._loops is None:
            self._loops = [
                e for e in range(self._element_length)
                if all(element[e] == 0 for element in self.cocircuits())
            ]
        return self._loops

    def plot(self, vertex_size: int = 600, figsize: int = None, aspect_ratio=None) -> None:
        r"""
        Plot the big face lattice of the oriented matroid.

        INPUT:

        For arguments, see :func:`.functions.plot_sign_vectors`.

        .. NOTE::

            Only works well for small length and dimension.
        """
        plot_sign_vectors(set().union(*self.faces_all()), vertex_size=vertex_size, figsize=figsize, aspect_ratio=aspect_ratio)

    @staticmethod
    def faces_from_vertices(vertices: set[SignVector], element_length: int) -> set[SignVector]:
        r"""
        Compute the covectors from the cocircuits.

        OUTPUT:

        - a set of all covectors of the oriented matroid.

        ALGORITHM:

        This function is based on an algorithm in [Fin01]_.

        .. [Fin01] Finschi, L.:
        „A graph theoretical approach for reconstruction and generation of oriented matroids“.
        PhD thesis. Zurich: ETH Zurich, 2001. doi: 10.3929/ethz-a-004255224.
        """
        covectors = {zero_sign_vector(element_length)}
        covectors_new = {zero_sign_vector(element_length)}
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

    @staticmethod
    def from_chirotopes(chirotopes: list[int] | str, rank: int, element_length: int) -> "OrientedMatroid":
        r"""
        Create an oriented matroid from chirotopes.

        INPUT:

        - ``chirotopes`` -- a list of integers or a string of ``+``, ``-``, ``0``.

        OUTPUT:

        - An instance of :class:`OrientedMatroid`.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: om = OrientedMatroid.from_chirotopes([1, 2, 3, 4, 6, 0], 2, 4)
            sage: om
            Oriented matroid of dimension 1 with elements of size 4.
            sage: om.faces_all()
            [{(0000)},
             {(+0--), (-0++), (--00), (++00), (0+++), (0---)},
             {(-+++), (+---), (--++), (++++), (----), (++--)}]

        ::

            sage: om = OrientedMatroid.from_chirotopes([0, 0, -1, -1, 0, 1, 1, 1, 1, 1], 2, 5)
            sage: om
            Oriented matroid of dimension 1 with elements of size 5.
            sage: om.faces_all()
            [{(00000)},
             {(-+++0), (+---0), (000++), (-++0-), (+--0+), (000--)},
             {(-+++-), (-++++), (+--++), (+----), (-++--), (+---+)}]

        Chirotopes can also be specified as a string::

            sage: OrientedMatroid.from_chirotopes("00--0+++++", 2, 5)
            Oriented matroid of dimension 1 with elements of size 5.
        """
        om = OrientedMatroid(rank=rank, element_length=element_length)
        for (indices, value) in zip(Combinations(element_length, rank), chirotopes):
            om.set_chirotope(indices, value)
        return om


class OrientedMatroidWithLattice(OrientedMatroid):
    r"""
    An oriented matroid that keeps track of the underlying lattice.
    """
    def __init__(self, matrix=None, rank: int = None, element_length: int = None) -> None:
        super().__init__(matrix=matrix, rank=rank, element_length=element_length)
        self.above = {}
        self.below = {}
        self._connected_dimensions = set()  # faces of this dimensions are already connected with lower faces

    def _connect(self, lower_face: SignVector, upper_face: SignVector):
        r"""Connect two faces."""
        if lower_face not in self.above:
            self.above[lower_face] = set()
        self.above[lower_face].add(upper_face)
        if upper_face not in self.below:
            self.below[upper_face] = set()
        self.below[upper_face].add(lower_face)

    def cocircuits(self):
        result = super().cocircuits()
        if 0 not in self._connected_dimensions:
            zero = zero_sign_vector(self._element_length)
            for face in result:
                self._connect(zero, face)
            self._connected_dimensions.add(0)
        return result

    def _lower_faces(self, dimension: int):
        if not self._faces_by_dimension.get(dimension):
            raise ValueError(f"Dimension {dimension} is not available. Available dimensions: {sorted(self._faces_by_dimension.keys())}.")
        if dimension == -1:
            return set()
        output = set()
        for same_support_faces in classes_same_support(self._faces_by_dimension[dimension]):
            p_classes = parallel_classes(same_support_faces, self._element_length)
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
                        if dimension not in self._connected_dimensions:
                            self._connect(lower_face, face)
                            self._connect(-lower_face, -face)
                            self._connect(lower_face, flipped_face)
                            self._connect(-lower_face, -flipped_face)
                same_support_faces.remove(-face)
        self._connected_dimensions.add(dimension)
        return output

    def _connect_missing(self):
        self.faces_all()
        for dimension in range(self.rank):
            if dimension not in self._connected_dimensions:
                self._lower_faces(dimension)

    def plot(self, vertex_size: int = 600, figsize: int = None, aspect_ratio=None) -> None:
        self._connect_missing()
        Poset(self.above).plot(
            vertex_size=vertex_size,
            element_color="white",
            vertex_shape="",
        ).show(figsize=figsize, aspect_ratio=aspect_ratio)

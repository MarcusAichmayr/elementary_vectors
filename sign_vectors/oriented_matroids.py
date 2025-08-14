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
from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors
from sign_vectors import sign_symbolic, SignVector, sign_vector, zero_sign_vector

from .utility import loops, classes_same_support, parallel_classes
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
        v = sign_symbolic(value)
        if v > 0:
            return cls.POS
        if v < 0:
            return cls.NEG
        return cls.ZERO


class OrientedMatroid(SageObject):
    r"""
    Class representing a (realizable) oriented matroid.

    EXAMPLES:

    We define an oriented matroid from a matrix::

        sage: from sign_vectors.oriented_matroids import *
        sage: M = matrix([[1, 2, 0, 0, 0], [0, 1, 2, 0, 0], [0, 0, 0, 1, 1]])
        sage: om = OrientedMatroid(M)
        sage: om
        OrientedMatroid of dimension 2 with covectors of size 5.
        sage: om.ground_set
        {0, 1, 2, 3, 4}
        sage: om.dimension
        2
        sage: om.rank
        3

    We can easily compute the cocircuits and topes::

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
         (000--),
         (--+--),
         (+0-++),
         (+----),
         (++000),
         (++0++),
         (---00),
         (++0--),
         (---++),
         (-0+00),
         (0++00),
         (-++00),
         (+--++),
         (+++00),
         (+++++),
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

    We can easily compute the faces of the same dimension::

        sage: om.faces(-1)
        {(00000)}
        sage: om.faces(0)
        {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)}
        sage: om.faces(1)
        {(+0---),
         (---00),
         (-0+++),
         (0----),
         (++0--),
         (++-00),
         (--+00),
         (-++00),
         (0++--),
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
        sage: om.all_faces()
        [{(00000)},
         {(--000), (000++), (+0-00), (0--00), (++000), (0++00), (000--), (-0+00)},
         {(+0---),
          (---00),
          (-0+++),
          (0----),
          (++0--),
          (++-00),
          (--+00),
          (-++00),
          (0++--),
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
         (---00),
         (-0+++),
         (0----),
         (++0--),
         (++-00),
         (--+00),
         (-++00),
         (0++--),
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
        sage: om.chirotopes_as_string()
        '0++++0++00'
        sage: om.chirotope([0, 1, 2])
        0
        sage: om.cocircuit([0, 1])
        (000++)
        sage: om.circuit([0, 1, 2, 3])
        (+-+00)

    We compute the dual oriented matroid::

        sage: om_dual = om.dual()
        sage: om_dual
        OrientedMatroid of dimension 1 with covectors of size 5.
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
        sage: om.all_faces()
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

        self._debug = False

    def _repr_(self) -> str:
        return f"OrientedMatroid of dimension {self.dimension} with covectors of size {self._element_length}."

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

        - The chirotope value as a Sign.

        .. SEEALSO::

            - :meth:`chirotopes`
            - :class:`Sign`

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
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
        self._chirotope_dict[tuple(indices)] = value

    def chirotopes(self) -> list[Sign]:
        r"""
        Compute all chirotopes of the oriented matroid.

        OUTPUT:

        - A list of chirotopes as Sign values.

        .. SEEALSO::

            - :meth:`chirotope`
            - :class:`Sign`

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotopes()
            [+, +, +, +, +, 0]
        """
        return [self.chirotope(indices) for indices in Combinations(self._element_length, self.rank)]

    def chirotopes_as_string(self) -> str:
        r"""Represent the chirotopes as a string."""
        return "".join(str(chirotope) for chirotope in self.chirotopes())

    def cocircuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute the cocircuit for the given indices.

        Cocircuits correspond to the elements in the row space of the matrix.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The cocircuit as a SignVector.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`circuit`

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
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
        Compute the circuit for the given indices.

        Circuits correspond to the elements in the kernel of the matrix.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The circuit as a SignVector.

        .. SEEALSO::

            - :meth:`circuits`
            - :meth:`cocircuit`

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
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

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.cocircuits()
            {(+0--), (-0++), (--00), (++00), (0+++), (0---)}
        """
        if 0 not in self._faces_by_dimension:
            if self._debug:
                print("Computing cocircuits...")
            result = set(self._cocircuit_generator())
            if result == set():
                return
            self._faces_by_dimension[0] = result
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

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.circuits()
            {(00-+), (+-0+), (-+0-), (00+-), (-+-0), (+-+0)}
        """
        if self._debug:
            print("Computing circuits...")
        return set(self._circuit_generator())

    def covectors(self) -> list[SignVector]:
        r"""
        Compute the covectors of the oriented matroid.

        .. NOTE::

            The covectors are all sign vectors corresponding to the row space.

        .. SEEALSO::

            - :meth:`cocircuits`
            - :meth:`topes`
            - :meth:`faces`
        """
        return self._covectors_from_cocircuits(self.cocircuits())

    def vectors(self) -> list[SignVector]:
        r"""
        Compute the vectors of the oriented matroid.

        .. NOTE::

            The vectors are all sign vectors corresponding to the kernel.

        .. SEEALSO::

            - :meth:`circuits`
        """
        return self._covectors_from_cocircuits(self.circuits())

    def topes(self) -> list[SignVector]:
        r"""
        Compute the topes of the oriented matroid.

        The topes are the support-maximal sign vectors corresponding to the row space.

        .. SEEALSO::

            - :meth:`cocircuits`
        """
        if self._faces_by_dimension.get(self.dimension) is None:
            if self._debug:
                print("Computing topes...")
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

            The results are hashed.

        .. SEEALSO::

            - :meth:`covectors`
            - :meth:`topes`
        """
        if self._debug:
            print(f"Faces available for: {self._faces_by_dimension.keys()}")
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

    def all_faces(self) -> list[set[SignVector]]:
        r"""
        Return all faces of the oriented matroid separated by dimension.

        .. SEEALSO::

            - :meth:`faces`
            - :func:`~face_enumeration`
        """
        return [self.faces(d) for d in range(-1, self.dimension + 1)]

    def _compute_lower_faces(self, dimension: int) -> None:
        if dimension - 1 not in self._faces_by_dimension:
            if self._debug:
                print(f"Computing faces for dimension {dimension - 1}...")
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
            - :func:`~face_enumeration`
        """
        if not self._faces_by_dimension.get(dimension):
            raise ValueError(f"Dimension {dimension} is not available. Available dimensions: {self._faces_by_dimension.keys()}.")
        output = set()
        for same_support_faces in classes_same_support(self._faces_by_dimension[dimension]):
            p_classes = parallel_classes(same_support_faces, self._element_length)
            for face in same_support_faces:
                for parallel_class in p_classes:
                    if all(face[i] == 0 for i in parallel_class):
                        continue
                    if face.reverse_signs_in(parallel_class) in same_support_faces:
                        output.add(sign_vector(0 if i in parallel_class else face[i] for i in range(self._element_length)))
        return output

    def _covectors_from_cocircuits(self, cocircuits: set[SignVector]) -> set[SignVector]:
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
        covectors = {zero_sign_vector(self._element_length)}
        covectors_new = {zero_sign_vector(self._element_length)}
        while covectors_new:
            element1 = covectors_new.pop()
            for element2 in cocircuits:
                if element2 <= element1:
                    continue
                new_element = element2.compose(element1)
                if new_element not in covectors:
                    covectors.add(new_element)
                    covectors_new.add(new_element)
        return covectors

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
            if self._debug:
                print("Computing loops...")
            self._loops = [
                e for e in range(self._element_length)
                if all(element[e] == 0 for element in self.cocircuits())
            ]
        return self._loops

    def plot(self, **kwargs) -> None:
        r"""
        Plot the big face lattice of the oriented matroid.

        .. NOTE::

            Only works well for small length and dimension.
        """
        plot_sign_vectors(set().union(*self.all_faces()), **kwargs)

    def dual(self) -> "OrientedMatroid":
        r"""Return the dual oriented matroid."""
        self.chirotopes() # compute all chirotopes
        om = OrientedMatroid(rank=self._element_length - self.rank, element_length=self._element_length)
        for indices, value in self._chirotope_dict.items():
            indices_set = set(indices)
            complement = tuple(i for i in range(self._element_length) if i not in indices_set)
            inversions = sum(i < j for i in indices for j in complement)
            om.set_chirotope(complement, -value if inversions & 1 else value) # check last bit
        return om


def cocircuits_from_matrix(M) -> set[SignVector]:
    r"""
    Compute a set of cocircuits determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix.

    OUTPUT:

    A set of cocircuits determined by the matrix ``M``.

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[1, 0, 2], [0, 1, -1]])
        sage: cocircuits_from_matrix(A)
        {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}
    """
    return OrientedMatroid(M).cocircuits()

def circuits_from_matrix(M) -> set[SignVector]:
    r"""
    Compute a set of circuits determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix.

    OUTPUT:

    A set of circuits determined by the matrix ``M``.

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[1, 0, 2], [0, 1, -1]])
        sage: circuits_from_matrix(A)
        {(-++), (+--)}
    """
    return OrientedMatroid(M).circuits()


def covectors_from_cocircuits(cocircuits: set[SignVector]) -> set[SignVector]:
    r"""
    Use an iterable of cocircuits to compute all covectors of the corresponding oriented matroid.

    INPUT:

    - ``cocircuits`` -- an iterable of cocircuits of an oriented matroid.

    OUTPUT:

    - a set of all covectors of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [Fin01]_.

    .. [Fin01] Finschi, L.:
       „A graph theoretical approach for reconstruction and generation of oriented matroids“.
       PhD thesis. Zurich: ETH Zurich, 2001. doi: 10.3929/ethz-a-004255224.

    EXAMPLES:

    First, we need cocircuits.
    For this purpose, we compute the cocircuits corresponding to some matrix::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[1, 2, 0], [0, 1, -1]])
        sage: ccA = cocircuits_from_matrix(A)
        sage: ccA
        {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}
        sage: covectors_from_cocircuits(ccA)
        {(000),
         (---),
         (0+-),
         (++-),
         (--0),
         (0-+),
         (--+),
         (++0),
         (+0+),
         (+++),
         (-0-),
         (-+-),
         (+-+)}
    """
    if not cocircuits:
        raise ValueError("List of cocircuits is empty.")
    for _ in cocircuits:
        length = _.length()
        break
    covectors = {zero_sign_vector(length)}
    covectors_new = {zero_sign_vector(length)}
    while covectors_new:
        element1 = covectors_new.pop()
        for element2 in cocircuits:
            if element2 <= element1:
                continue
            new_element = element2.compose(element1)
            if new_element not in covectors:
                covectors.add(new_element)
                covectors_new.add(new_element)
    return covectors


def topes_from_cocircuits(cocircuits: set[SignVector]) -> set[SignVector]:
    r"""
    Use the cocircuits of an oriented matroid to compute the topes.

    INPUT:

    - ``cocircuits`` -- an iterable of cocircuits of an oriented matroid.

    OUTPUT:

    A set of topes of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [Fin01]_.

    EXAMPLES:

    First, we need cocircuits.
    For this purpose, we compute the cocircuits corresponding to some matrix::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[1, 0, 2], [0, 1, -1]])
        sage: ccA = cocircuits_from_matrix(A)
        sage: ccA
        {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}
        sage: topes_from_cocircuits(ccA)
        {(---), (++-), (--+), (+++), (-+-), (+-+)}
    """
    if not cocircuits:
        raise ValueError("List is empty.")
    for _ in cocircuits:
        length = _.length()
        break

    covectors = {zero_sign_vector(length)}
    covectors_new = {zero_sign_vector(length)}
    topes = set()
    loop_list = loops(cocircuits)

    while covectors_new:
        element1 = covectors_new.pop()
        for element2 in cocircuits:
            if element2 <= element1:
                continue
            new_element = element2.compose(element1)
            if new_element not in covectors:
                covectors.add(new_element)
                if new_element.zero_support() == loop_list:
                    topes.add(new_element)
                else:
                    covectors_new.add(new_element)
    return topes


def lower_faces(i_faces: set[SignVector]) -> set[SignVector]:
    r"""
    Compute the faces of lower dimension.

    INPUT:

    - ``i_faces`` -- an iterable of all faces with the same dimension ``i`` of an oriented matroid.

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

        - :func:`~face_enumeration`
        - :func:`~covectors_from_topes`
    """
    if not i_faces:
        raise ValueError("List is empty.")
    for _ in i_faces:
        length = _.length()
        break
    output = set()
    for same_support_faces in classes_same_support(i_faces):
        p_classes = parallel_classes(same_support_faces, length)
        for face in same_support_faces:
            for parallel_class in p_classes:
                if all(face[i] == 0 for i in parallel_class):
                    continue
                if face.reverse_signs_in(parallel_class) in same_support_faces:
                    output.add(sign_vector(0 if i in parallel_class else face[i] for i in range(length)))
    return output


def face_enumeration(i_faces: set[SignVector]) -> list[set[SignVector]]:
    r"""
    Compute all faces with smaller dimension.

    INPUT:

    - ``i_faces`` -- an iterable of all faces with the same dimension ``i`` of an oriented matroid.

    OUTPUT:
    A list of sets. Every set consists of all faces of the same dimension
    smaller than or equal to ``i`` of the oriented matroid.

    Apply this function to the topes and obtain all faces of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [FST91]_.
    See also [Fin01]_.

    .. SEEALSO::

        - :func:`~lower_faces`
        - :func:`~covectors_from_topes`
        - :func:`~covectors_from_matrix`

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[1, 0, 2], [0, 1, -1]])
        sage: topes = topes_from_cocircuits(cocircuits_from_matrix(A))
        sage: topes
        {(---), (++-), (--+), (+++), (-+-), (+-+)}
        sage: face_enumeration(topes)
        [{(000)},
         {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
         {(---), (-+-), (++-), (--+), (+++), (+-+)}]
    """
    if not i_faces:
        raise ValueError("List is empty.")
    faces = [set(i_faces)]

    while len(faces[0]) > 1:
        faces.insert(0, lower_faces(faces[0]))
    return faces

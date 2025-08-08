r"""
Oriented matroids
=================

We define some matrix::

    sage: A = matrix([[2, -1, -1]])
    sage: A
    [ 2 -1 -1]

Cocircuits
~~~~~~~~~~

Now, we compute the cocircuits of the oriented matroid corresponding to the kernel
of the matrix ``A``::

    sage: from sign_vectors.oriented_matroids import *
    sage: ccA = cocircuits_from_matrix(A)
    sage: ccA
    {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}

Covectors
~~~~~~~~~

We can also use the cocircuits to compute all covectors of the corresponding
oriented matroid::

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
    sage: covectors_from_matrix(A)
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

Topes
~~~~~

Next, we compute the topes::

    sage: tA = topes_from_cocircuits(ccA)
    sage: tA
    {(---), (++-), (--+), (+++), (-+-), (+-+)}
    sage: topes_from_matrix(A)
    {(---), (++-), (--+), (+++), (-+-), (+-+)}

Face enumeration
~~~~~~~~~~~~~~~~

Next, we compute all covectors separated by their rank::

    sage: face_enumeration(tA)
    [{(000)},
     {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
     {(---), (-+-), (++-), (--+), (+++), (+-+)}]

By passing ``dual=False``, we can compute the covectors of the
dual oriented matroid::

    sage: cocircuits_from_matrix(A, dual=False)
    {(+--), (-++)}
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

from collections.abc import Generator
from enum import IntEnum

from sage.combinat.combination import Combinations
from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors
from sign_vectors import sign_symbolic, SignVector, sign_vector, zero_sign_vector

from .utility import loops, classes_same_support, parallel_classes


class Sign(IntEnum):
    r"""
    Sign values used in oriented matroids.

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
        """Return the opposite sign."""
        return Sign(-self.value) if self.value != 0 else Sign.ZERO

    @classmethod
    def _missing_(cls, value):
        v = sign_symbolic(value)
        if v > 0:
            return cls.POS
        if v < 0:
            return cls.NEG
        return cls.ZERO


# TODO duplicates code of ElementaryVectors
class OrientedMatroid(SageObject):
    r"""
    Class representing the chirotopes of an oriented matroid.
    """
    def __init__(self, M, hash_faces: bool = True) -> None:
        try:
            self.matrix = M.matrix_from_rows(M.pivot_rows())
        except NotImplementedError as exc:
            if all(minor == 0 for minor in M.minors(M.nrows())):
                raise ValueError("Provide a matrix with maximal rank.") from exc
            self.matrix = M
        self.rank, self.length = self.matrix.dimensions()
        self._chirotopes = {}
        self._hash_faces = hash_faces
        self._faces = {0: set(zero_sign_vector(self.length))}
        self._faces_dual = {0: set(zero_sign_vector(self.length))}
        self._loops = None
        self._coloops = None

    def chirotope(self, indices: list[int]) -> Sign:
        r"""
        Compute the chirotope for the given indices.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The chirotope value as a Sign.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotope([1, 2])
            +
        """
        indices = tuple(indices)
        chirotope = self._chirotopes.get(indices)
        if chirotope is None:
            try:
                chirotope = Sign(self.matrix.matrix_from_columns(indices).det())
                self._chirotopes[indices] = chirotope
            except ValueError as e:
                raise ValueError(f"Indices {indices} should have size {self.rank} and not {len(indices)}.") from e
        return chirotope

    def chirotopes(self) -> list[Sign]:
        r"""
        Compute all chirotopes of the oriented matroid.

        OUTPUT:

        - A list of chirotopes as Sign values.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.chirotopes()
            [+, +, +, +, +, 0]
        """
        return [self.chirotope(indices) for indices in Combinations(self.length, self.rank)]

    def circuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute the circuit for the given indices.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The circuit as a SignVector.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.circuit([0, 1, 2])
            (+-+0)
        """
        element = [0] * self.length
        for pos in range(self.rank + 1):
            indices_chirotope = indices.copy()
            i = indices_chirotope.pop(pos)
            chirotope = self.chirotope(indices_chirotope)
            if chirotope != 0:
                # check oddness of last bit of pos
                element[i] = -chirotope if (pos & 1) else chirotope
        return sign_vector(element)

    def cocircuit(self, indices: list[int]) -> SignVector:
        r"""
        Compute the cocircuit for the given indices.

        INPUT:

        - ``indices`` -- a list of indices.

        OUTPUT:

        - The cocircuit as a SignVector.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.cocircuit([0])
            (0---)
        """
        element = [0] * self.length
        pos = 0
        for i in range(self.length):
            if i in indices:
                pos += 1
                continue
            chirotope = self.chirotope(sorted(indices + [i]))
            if chirotope != 0:
                # check oddness of last bit of pos
                element[i] = -chirotope if (pos & 1) else chirotope
        return sign_vector(element)

    def circuits_generator(self) -> Generator[SignVector]:
        r"""
        Compute the circuits of the oriented matroid.

        OUTPUT:

        - A generator of circuits as SignVectors.
        """
        for indices in Combinations(self.length, self.rank + 1):
            yield self.circuit(indices)
            yield -self.circuit(indices)

    def circuits(self) -> set[SignVector]:
        r"""
        Compute the circuits of the oriented matroid.

        OUTPUT:

        - A set of circuits as SignVectors.

        .. NOTE::

            The circuits are the support-minimal sign vectors corresponding to the kernel.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.circuits()
            {(00-+), (+-0+), (-+0-), (00+-), (-+-0), (+-+0)}
        """
        if 1 in self._faces:
            return self._faces[1]
        circuits = set(self.circuits_generator())
        if self._hash_faces:
            self._faces[1] = circuits
        return circuits

    def cocircuits_generator(self) -> Generator[SignVector]:
        r"""
        Compute the cocircuits of the oriented matroid.

        OUTPUT:

        - A generator of cocircuits as SignVectors.
        """
        for indices in Combinations(self.length, self.rank - 1):
            yield self.cocircuit(indices)
            yield -self.cocircuit(indices)

    def cocircuits(self) -> set[SignVector]:
        r"""
        Compute the cocircuits of the oriented matroid.

        OUTPUT:

        - A set of cocircuits as SignVectors.

        .. NOTE::

            The cocircuits are the support-minimal sign vectors corresponding to the row space.

        EXAMPLES::

            sage: from sign_vectors.oriented_matroids import *
            sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
            sage: om = OrientedMatroid(M)
            sage: om.cocircuits()
            {(+0--), (-0++), (--00), (++00), (0+++), (0---)}
        """
        if 1 in self._faces_dual:
            return self._faces_dual[1]
        cocircuits = set(self.cocircuits_generator())
        if self._hash_faces:
            self._faces_dual[1] = cocircuits
        return cocircuits

    def vectors(self) -> list[SignVector]:
        r"""
        Compute the vectors of the oriented matroid.

        .. NOTE::

            The vectors are all sign vectors corresponding to the kernel.
        """
        # TODO should store support-maximal elements
        return covectors_from_cocircuits(self.circuits())

    def covectors(self) -> list[SignVector]:
        r"""
        Compute the covectors of the oriented matroid.

        .. NOTE::

            The covectors are all sign vectors corresponding to the row space.
        """
        # TODO should store topes
        return covectors_from_cocircuits(self.cocircuits())

    def support_maximal_elements(self) -> list[SignVector]:
        r"""
        Compute the support-maximal elements of the oriented matroid.

        OUTPUT:

        - A list of support-maximal elements as SignVectors.
        """
        return topes_from_cocircuits(self.circuits())

    def topes(self) -> list[SignVector]:
        r"""
        Compute the topes of the oriented matroid.

        .. NOTE::

            The topes are the support-maximal sign vectors corresponding to the row space.
        """
        return topes_from_cocircuits(self.cocircuits())

    def loops(self) -> list[int]:
        r"""
        Compute the loops of the oriented matroid.

        Output:
        A list of zero entries of each covector.
        """
        if self._loops is not None:
            return self._loops
        # TODO we check twice as many elements as necessary
        self._loops = [
            i for i in range(self.length)
            if all(circuit[i] == 0 for circuit in self.circuits())
        ]
        return self._loops

    def coloops(self) -> list[int]:
        r"""
        Compute the coloops of the oriented matroid.

        Output:
        A list of zero entries of each circuit.
        """
        if self._coloops is not None:
            return self._coloops
        self._coloops = [i for i, column in enumerate(self.matrix.columns()) if column == 0]
        return self._coloops

    def faces(self, level: int):
        r"""
        Compute the faces of the same level of the oriented matroid.
        """
        if level in self._faces:
            return self._faces[level]
        if level == 1:
            return self.circuits()
        if level == self.rank:
            return self.support_maximal_elements()
        # should use lower faces to get down
        raise NotImplementedError

    def dual_faces(self, level: int):
        r"""
        Compute the dual faces of the same level of the oriented matroid.
        """
        raise NotImplementedError

    def dimension(self):
        raise NotImplementedError

    def plot(self):
        r"""Plot the big face lattice of the oriented matroid."""
        raise NotImplementedError


# TODO redundant
class Cocircuits(ElementaryVectors):
    r"""
    Class used to compute cocircuits and circuits.

    TESTS::

        sage: from sign_vectors.oriented_matroids import *
        sage: M = matrix([[1, 2, 4, 0], [0, 1, 2, 0]])
        sage: cc = Cocircuits(M)
        sage: cc.element([1, 2, 3])
        Traceback (most recent call last):
        ...
        ValueError: The indices [1, 2, 3] correspond to the zero vector!
        sage: cc.elements()
        {(0+-0), (000-), (0-+0), (000+)}
        sage: cc.minors
        {(0, 1): 1, (0, 2): 1, (0, 3): 0, (1, 2): 0, (1, 3): 0, (2, 3): 0}
    """
    def _compute_minor(self, indices: tuple[int]) -> int:
        return sign_symbolic(super()._compute_minor(indices))

    def _zero_element(self) -> list[int]:
        return [0] * self.length

    def _element_kernel(self, indices: list[int], mark_zeros: bool = False) -> SignVector:
        return sign_vector(super()._element_kernel(indices, mark_zeros=mark_zeros))

    def _element_row_space(self, indices: list[int], mark_zeros: bool = False) -> SignVector:
        return sign_vector(super()._element_row_space(indices, mark_zeros=mark_zeros))

    def generator(
        self,
        dual: bool = True,
        prevent_multiples: bool = True,
        reverse: bool = False
    ) -> Generator[SignVector]:
        for cocircuit in super().generator(dual=dual, prevent_multiples=prevent_multiples, reverse=reverse):
            yield cocircuit
            yield -cocircuit

    def elements(self, dual: bool = True, prevent_multiples: bool = True) -> set[SignVector]:
        return set(self.generator(dual=dual, prevent_multiples=prevent_multiples))


def cocircuits_from_matrix(M, dual: bool = True):
    r"""
    Compute a set of cocircuits determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix with real arguments.

    - ``dual`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``dual`` is true, returns a set of cocircuits determined by the
      kernel of the matrix ``M``.

    - If ``dual`` is false, returns a set of cocircuits determined by the row
      space of the matrix ``M``.

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: cocircuits_from_matrix(A)
        {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}
        sage: B = matrix([[1, 0, 0, 0], [0, 1, 0, 2], [0, 0, 1, -1]])
        sage: B
        [ 1  0  0  0]
        [ 0  1  0  2]
        [ 0  0  1 -1]
        sage: cocircuits_from_matrix(B, dual=False)
        {(-000), (00-+), (0-0-), (0--0), (0+0+), (00+-), (0++0), (+000)}
    """
    return Cocircuits(M).elements(dual=dual)


def covectors_from_cocircuits(cocircuits):
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
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
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
        covector1 = covectors_new.pop()
        for covector2 in cocircuits:
            if covector2 <= covector1:
                continue
            composition = covector2.compose(covector1)
            if composition not in covectors:
                covectors.add(composition)
                covectors_new.add(composition)
    return covectors


def topes_from_cocircuits(cocircuits):
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
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
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
        covector1 = covectors_new.pop()
        for covector2 in cocircuits:
            if covector2 <= covector1:
                continue
            composition = covector2.compose(covector1)
            if composition not in covectors:
                covectors.add(composition)
                if composition.zero_support() == loop_list:
                    topes.add(composition)
                else:
                    covectors_new.add(composition)
    return topes


def lower_faces(covectors):
    r"""
    Compute the lower faces of given covectors.

    INPUT:

    - ``covectors`` -- an iterable of all covectors with same rank ``r`` of an oriented matroid.

    OUTPUT:

    Returns a set of covectors of rank ``r-1`` of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [FST91]_.
    See also [Fin01]_.

    .. [FST91] Fukuda, K., Saito, S., and Tamura, A.:
       "Combinatorial face enumeration in arrangements and oriented matroids".
       In: Discrete Applied Mathematics 31.2 (1991), pp. 141-149.
       doi: 10.1016/0166-218X(91)90066-6.

    .. SEEALSO::

        :func:`~face_enumeration`
        :func:`~covectors_from_topes`
    """
    if not covectors:
        raise ValueError("List is empty.")
    for _ in covectors:
        length = _.length()
        break
    output = set()
    for covectors_with_same_support in classes_same_support(covectors):
        p_classes = parallel_classes(covectors_with_same_support)
        for covector in covectors_with_same_support:
            for parallel_class in p_classes:
                for i in parallel_class:
                    if not covector[i]:
                        continue
                    if (
                        covector.reverse_signs_in(parallel_class)
                        in covectors_with_same_support
                    ):
                        output.add(
                            sign_vector(
                                0 if i in parallel_class else covector[i]
                                for i in range(length)
                            )
                        )
                    break
    return output


def face_enumeration(covectors):
    r"""
    Compute all covectors with less rank than the given covectors.

    INPUT:

    - ``covectors`` -- an iterable of all covectors of same rank ``r`` of an oriented matroid.

    OUTPUT:

    Returns a list of sets. Every set consists of all covectors of the same rank
    smaller than or equal to ``r`` of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [FST91]_.
    See also [Fin01]_.

    .. SEEALSO::

        :func:`~lower_faces`
        :func:`~covectors_from_topes`
        :func:`~covectors_from_matrix`

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        {(---), (++-), (--+), (+++), (-+-), (+-+)}
        sage: face_enumeration(tA)
        [{(000)},
         {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
         {(---), (-+-), (++-), (--+), (+++), (+-+)}]
    """
    if not covectors:
        raise ValueError("List is empty.")
    faces = [set(covectors)]

    while len(faces[0]) > 1:
        faces.insert(0, lower_faces(faces[0]))
    return faces


def topes_from_matrix(M, dual: bool = True):
    r"""
    Return the topes of the oriented matroid corresponding to the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``dual`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``dual`` is true, returns a set of topes determined by the kernel of
      the matrix ``M``.

    - If ``dual`` is false, returns a set of topes determined by the row space
      of the matrix ``M``.

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: topes_from_matrix(A)
        {(---), (++-), (--+), (+++), (-+-), (+-+)}

    TESTS::

        sage: M = zero_matrix(1, 3)
        sage: topes_from_matrix(M, dual=False)
        {(000)}
    """
    cocircuits = cocircuits_from_matrix(M, dual=dual)
    if cocircuits:
        return topes_from_cocircuits(cocircuits)
    return {zero_sign_vector(M.ncols())}


def covectors_from_topes(topes, separate: bool = False):
    r"""
    Compute all covectors from the topes.

    INPUT:

    - ``topes`` -- an iterable of topes.

    - ``separate`` -- a boolean (default: ``False``)

    OUTPUT:

    The set of covectors of the corresponding oriented matroid.

    - If ``separate`` is false, returns a list of covectors. The covectors are
      sorted by rank. (default)

    - If ``separate`` is true, returns a list of lists of covectors, separated
      by their rank.

    .. SEEALSO::

        :func:`~face_enumeration`

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        {(---), (++-), (--+), (+++), (-+-), (+-+)}
        sage: covectors_from_topes(tA)
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
        sage: covectors_from_topes(tA, separate=True)
        [{(000)},
         {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
         {(---), (-+-), (++-), (--+), (+++), (+-+)}]
    """
    if separate:
        return face_enumeration(topes)
    return set().union(*face_enumeration(topes))


def cocircuits_from_topes(topes):
    r"""
    Compute all cocircuits from the topes.

    INPUT:

    - ``topes`` -- an iterable of topes.

    OUTPUT:

    A set of cocircuits of the corresponding oriented matroid.

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import *
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        {(---), (++-), (--+), (+++), (-+-), (+-+)}
        sage: cocircuits_from_topes(tA)
        {(0+-), (--0), (0-+), (++0), (+0+), (-0-)}
    """
    return face_enumeration(topes)[1]


def covectors_from_matrix(M, dual: bool = True, algorithm: str = None, separate: bool = False):
    r"""
    Return the covectors of the oriented matroid corresponding to the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix.

    - ``dual`` -- a boolean (default: ``True``)

    - ``algorithm`` -- either ``None`` (default), ``"face_enumeration"`` or ``"fe"``

    - ``separate`` -- a boolean (default: ``False``)

    OUTPUT:

    Returns a set of covectors of an oriented matroid corresponding to the
    matrix ``M``.

    - If ``dual`` is true, the returned covectors will be determined by the
      kernel of the matrix ``M``.

    - If ``dual`` is false, the returned covectors will be determined by the
      row space of the matrix ``M``.

    - If ``algorithm`` is ``"face_enumeration"`` or the shortcut ``"fe"``,
      applies the algorithm face enumeration.

    - If ``separate`` is true, returns a list of sets of covectors, separated
      by their rank by applying the algorithm face enumeration.

    - If ``separate`` is false, returns a set of covectors.

    .. SEEALSO::

        :func:`~face_enumeration`

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import *

    We define some matrix::

        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: covectors_from_matrix(A)
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
        sage: covectors_from_matrix(A, separate=True)
        [{(000)},
         {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
         {(---), (-+-), (++-), (--+), (+++), (+-+)}]
        sage: covectors_from_matrix(A, algorithm="face_enumeration", separate=True)
        [{(000)},
         {(0+-), (--0), (0-+), (++0), (+0+), (-0-)},
         {(---), (-+-), (++-), (--+), (+++), (+-+)}]

    TESTS::

        sage: covectors_from_matrix(zero_matrix(1, 4), dual=False, separate=False)
        {(0000)}
        sage: covectors_from_matrix(zero_matrix(1, 4), dual=False, separate=True)
        [{(0000)}]
    """
    if algorithm is None:
        if separate:
            algorithm = "face_enumeration"
        else:
            cocircuits = cocircuits_from_matrix(M, dual=dual)
            if cocircuits:
                return covectors_from_cocircuits(cocircuits)
            return {zero_sign_vector(M.ncols())}
    if algorithm in ["face_enumeration", "fe"]:
        return covectors_from_topes(
            topes_from_matrix(M, dual=dual), separate=separate
        )
    raise ValueError(f"no algorithm '{algorithm}'")

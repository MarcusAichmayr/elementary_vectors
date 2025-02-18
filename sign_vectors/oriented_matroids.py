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
of the matrix ``A``.
(Cocircuits are minimal nonzero elements of an oriented matroid
with respect to the conformal relation.)::

    sage: from sign_vectors.oriented_matroids import *
    sage: ccA = cocircuits_from_matrix(A)
    sage: ccA
    {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}

Covectors
~~~~~~~~~

We can also use the cocircuits to compute all covectors of the corresponding
oriented matroid::

    sage: covectors_from_cocircuits(ccA)
    {(000),
     (+-+),
     (---),
     (-+-),
     (0-+),
     (+0+),
     (--0),
     (-0-),
     (+++),
     (0+-),
     (++-),
     (--+),
     (++0)}

Topes
~~~~~

Next, we compute the topes using the cocircuits.
(Topes are the covectors that are maximal with respect to the conformal relation)::

    sage: tA = topes_from_cocircuits(ccA)
    sage: tA
    {(+-+), (---), (-+-), (+++), (--+), (++-)}

There are some further commands to work with oriented matroids::

    sage: covectors_from_matrix(A)
    {(000),
     (+-+),
     (---),
     (-+-),
     (0-+),
     (+0+),
     (--0),
     (-0-),
     (+++),
     (0+-),
     (++-),
     (--+),
     (++0)}
    sage: topes_from_matrix(A)
    {(+-+), (---), (-+-), (+++), (--+), (++-)}
    sage: covectors_from_topes(tA)
    {(000),
     (+-+),
     (---),
     (-+-),
     (0-+),
     (+0+),
     (--0),
     (-0-),
     (++-),
     (+++),
     (0+-),
     (--+),
     (++0)}
    sage: cocircuits_from_topes(tA)
    {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}

Face enumeration
~~~~~~~~~~~~~~~~

Next, we compute all covectors separated by their rank::

    sage: face_enumeration(tA)
    [{(000)},
     {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
     {(+-+), (---), (-+-), (--+), (++-), (+++)}]
    sage: covectors_from_matrix(A, algorithm="face_enumeration", separate=True)
    [{(000)},
     {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
     {(+-+), (---), (-+-), (--+), (++-), (+++)}]

By passing ``dual=False``, we can compute the covectors of the
dual oriented matroid::

    sage: cocircuits_from_matrix(A, dual=False)
    {(-++), (+--)}
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

from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors
from sign_vectors import Sign, SignVector, sign_vector, zero_sign_vector

from .utility import loops, classes_same_support, parallel_classes


class Cocircuits(ElementaryVectors):
    r"""Class used to compute cocircuits and circuits."""
    def minor(self, indices):
        indices = tuple(indices)
        try:
            return super().minor(indices)
        except KeyError:
            self.minors[indices] = Sign(self.matrix.matrix_from_columns(indices).det())
            return self.minors[indices]

    def _zero_element(self) -> list:
        return [0] * self.length

    def element_kernel(self, indices: list, prevent_multiple: bool = False) -> SignVector:
        return sign_vector(super().element_kernel(indices, prevent_multiple=prevent_multiple))

    def element_row_space(self, indices: list, prevent_multiple: bool = False) -> SignVector:
        return sign_vector(super().element_row_space(indices, prevent_multiple=prevent_multiple))

    def generator(
        self,
        dual: bool = True,
        prevent_multiples: bool = True,
        reverse: bool = False
    ) -> Generator:
        for cocircuit in super().generator(dual=dual, prevent_multiples=prevent_multiples, reverse=reverse):
            yield cocircuit
            yield -cocircuit

    def elements(self, dual: bool = True, prevent_multiples: bool = True) -> set:
        return set(self.generator(dual=dual, prevent_multiples=prevent_multiples))


class OrientedMatroid(SageObject):
    def __init__(self, M) -> None:
        self.matrix = M
        self._cocircuits = Cocircuits(M)

    def circuits(self):
        return self._cocircuits.elements(dual=True)

    def cocircuits(self):
        return self._cocircuits.elements(dual=False)

    def vectors(self):
        return covectors_from_cocircuits(self.circuits())

    def covectors(self):
        return covectors_from_cocircuits(self.cocircuits())

    def topes(self):
        return topes_from_cocircuits(self.circuits())

    def cotopes(self):
        # TODO name?
        return topes_from_cocircuits(self.cocircuits())

    def faces(self, level: int):
        # should use lower faces to get down
        raise NotImplementedError

    def rank(self):
        raise NotImplementedError

    def dimension(self):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

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
        {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}
        sage: B = matrix([[1, 0, 0, 0], [0, 1, 0, 2], [0, 0, 1, -1]])
        sage: B
        [ 1  0  0  0]
        [ 0  1  0  2]
        [ 0  0  1 -1]
        sage: cocircuits_from_matrix(B, dual=False)
        {(0+0+), (00+-), (-000), (00-+), (0-0-), (0--0), (+000), (0++0)}
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
        {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}
        sage: covectors_from_cocircuits(ccA)
        {(000),
         (+-+),
         (---),
         (-+-),
         (0-+),
         (+0+),
         (--0),
         (-0-),
         (+++),
         (0+-),
         (++-),
         (--+),
         (++0)}
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
        {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}
        sage: topes_from_cocircuits(ccA)
        {(+-+), (---), (-+-), (+++), (--+), (++-)}
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
       „Combinatorial face enumeration in arrangements and oriented matroids“.
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
        {(+-+), (---), (-+-), (+++), (--+), (++-)}
        sage: face_enumeration(tA)
        [{(000)},
         {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
         {(+-+), (---), (-+-), (--+), (++-), (+++)}]
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
        {(+-+), (---), (-+-), (+++), (--+), (++-)}

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
        {(+-+), (---), (-+-), (+++), (--+), (++-)}
        sage: covectors_from_topes(tA)
        {(000),
         (+-+),
         (---),
         (-+-),
         (0-+),
         (+0+),
         (--0),
         (-0-),
         (++-),
         (+++),
         (0+-),
         (--+),
         (++0)}
        sage: covectors_from_topes(tA, separate=True)
        [{(000)},
         {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
         {(+-+), (---), (-+-), (--+), (++-), (+++)}]
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
        {(+-+), (---), (-+-), (+++), (--+), (++-)}
        sage: cocircuits_from_topes(tA)
        {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}
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
         (+-+),
         (---),
         (-+-),
         (0-+),
         (+0+),
         (--0),
         (-0-),
         (+++),
         (0+-),
         (++-),
         (--+),
         (++0)}
        sage: covectors_from_matrix(A, separate=True)
        [{(000)},
         {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
         {(+-+), (---), (-+-), (--+), (++-), (+++)}]
        sage: covectors_from_matrix(A, algorithm="face_enumeration", separate=True)
        [{(000)},
         {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
         {(+-+), (---), (-+-), (--+), (++-), (+++)}]
        sage: covectors_from_matrix(A, algorithm="fe", separate=True)
        [{(000)},
         {(0-+), (+0+), (--0), (-0-), (0+-), (++0)},
         {(+-+), (---), (-+-), (--+), (++-), (+++)}]

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

r"""
Computing with oriented matroids.

This module is about computation with oriented matroids.

EXAMPLES::

    sage: from sign_vectors.oriented_matroids import *

We define some matrix::

    sage: A = matrix([[1,2,0],[0,1,-1]])
    sage: A
    [ 1  2  0]
    [ 0  1 -1]

Now, we compute the cocircuits of the oriented matroid corresponding to the rows
of the matrix ``A``.
(Cocircuits are minimal non-zero elements of an oriented matroid
with respect to the conformal relation.)::

    sage: ccA = cocircuits_from_matrix(A, kernel=False)
    sage: ccA
    [(--0), (++0), (-0-), (+0+), (0-+), (0+-)]

We can also use the cocircuits to compute all covectors of the corresponding
oriented matroid::

    sage: covectors_from_cocircuits(ccA)
    [(---),
     (--+),
     (-0-),
     (000),
     (0-+),
     (--0),
     (+-+),
     (+0+),
     (+++),
     (++-),
     (-+-),
     (0+-),
     (++0)]

Next, we compute the topes using the cocircuits.
(Topes are the covectors that are maximal with respect to the conformal relation)::

    sage: tA = topes_from_cocircuits(ccA)
    sage: tA
    [(--+), (+++), (---), (+-+), (-+-), (++-)]

There are some further commands to work with oriented matroids::

    sage: covectors_from_matrix(A, kernel=False)
    [(---),
     (--+),
     (-0-),
     (000),
     (0-+),
     (--0),
     (+-+),
     (+0+),
     (+++),
     (++-),
     (-+-),
     (0+-),
     (++0)]
    sage: topes_from_matrix(A, kernel=False)
    [(--+), (+++), (---), (+-+), (-+-), (++-)]
    sage: covectors_from_topes(tA)
    [(000),
     (-0-),
     (0-+),
     (--0),
     (+0+),
     (0+-),
     (++0),
     (--+),
     (+++),
     (---),
     (+-+),
     (-+-),
     (++-)]
    sage: cocircuits_from_topes(tA)
    [(-0-), (0-+), (--0), (+0+), (0+-), (++0)]

By passing ``kernel=True`` (default), we can compute the covectors of the
dual oriented matroid::

    sage: cocircuits_from_matrix(A, kernel=True)
    [(-++), (+--)]
    sage: cocircuits_from_matrix(A)
    [(-++), (+--)]
    sage: topes_from_matrix(A, kernel=True)
    [(-++), (+--)]
    sage: covectors_from_matrix(A, kernel=True)
    [(+--), (-++), (000)]

Next, we compute all covectors separated by their rank::

    sage: face_enumeration(tA)
    [[(000)],
     [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
     [(--+), (+++), (---), (+-+), (-+-), (++-)]]
    sage: covectors_from_matrix(A, kernel=False, algorithm="face_enumeration", separate=True)
    [[(000)],
     [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
     [(--+), (+++), (---), (+-+), (-+-), (++-)]]
"""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.misc.flatten import flatten
from elementary_vectors import elementary_vectors
from sign_vectors import sign_vector, zero_sign_vector
from sign_vectors.utility import loops, classes_same_support, parallel_classes


def cocircuits_from_matrix(M, kernel=True):
    r"""
    Compute a list of cocircuits determined by the matrix ``A``.

    INPUT:

    - ``M`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``kernel`` is true, returns a list of cocircuits determined by the
      kernel of the matrix ``M``.

    - If ``kernel`` is false, returns a list of cocircuits determined by the row
      space of the matrix ``M``.

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: cocircuits_from_matrix(A)
        [(--0), (++0), (-0-), (+0+), (0-+), (0+-)]
        sage: B = matrix([[1,0,0,0],[0,1,0,2],[0,0,1,-1]])
        sage: B
        [ 1  0  0  0]
        [ 0  1  0  2]
        [ 0  0  1 -1]
        sage: cocircuits_from_matrix(B, kernel=False)
        [(+000), (-000), (0--0), (0++0), (0-0-), (0+0+), (00-+), (00+-)]
    """
    return sum([
        [sign_vector(v), sign_vector(-v)]
        for v in elementary_vectors(M, kernel=kernel)
    ], [])


def covectors_from_cocircuits(cocircuits):
    r"""
    Use a list of cocircuits to compute all covectors of the corresponding oriented matroid.

    INPUT:

    - ``cocircuits`` -- a list of cocircuits of an oriented matroid.

    OUTPUT:

    - a list of all covectors of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [Fin01]_.

    .. [Fin01] Finschi, L.:
       „A graph theoretical approach for reconstruction and generation of oriented matroids“.
       PhD thesis. Zurich: ETH Zurich, 2001. doi: 10.3929/ethz-a-004255224.

    EXAMPLES:

    First, we need a list of cocircuits.
    For this purpose, we compute the cocircuits corresponding to some matrix::

        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix, covectors_from_cocircuits
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: ccA = cocircuits_from_matrix(A)
        sage: ccA
        [(--0), (++0), (-0-), (+0+), (0-+), (0+-)]
        sage: covectors_from_cocircuits(ccA)
        [(---),
         (--+),
         (-0-),
         (000),
         (0-+),
         (--0),
         (+-+),
         (+0+),
         (+++),
         (++-),
         (-+-),
         (0+-),
         (++0)]
    """
    if not cocircuits:
        raise ValueError('List of cocircuits is empty.')
    n = cocircuits[0].length()
    covectors = {zero_sign_vector(n)}
    covectors_new = {zero_sign_vector(n)}
    while covectors_new != set():
        Y = covectors_new.pop()
        for X in cocircuits:
            if not X <= Y:  # otherwise Z = X.compose(Y) = Y in ``covectors``
                Z = X.compose(Y)
                if Z not in covectors:
                    covectors.add(Z)
                    covectors_new.add(Z)
    return list(covectors)


def topes_from_cocircuits(cocircuits):
    r"""
    Use the cocircuits of an oriented matroid to compute the topes.

    INPUT:

    - ``cocircuits`` -- a list of cocircuits of an oriented matroid.

    OUTPUT:

    A list of topes of the oriented matroid.

    ALGORITHM:

    This function is based on an algorithm in [Fin01]_.

    EXAMPLES:

    First, we need a list of cocircuits.
    For this purpose, we compute the cocircuits corresponding to some matrix::

        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix, topes_from_cocircuits
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: ccA = cocircuits_from_matrix(A)
        sage: ccA
        [(--0), (++0), (-0-), (+0+), (0-+), (0+-)]
        sage: topes_from_cocircuits(ccA)
        [(--+), (+++), (---), (+-+), (-+-), (++-)]
    """
    if not cocircuits:
        raise ValueError('List is empty.')
    n = cocircuits[0].length()

    covectors = {zero_sign_vector(n)}
    covectors_new = {zero_sign_vector(n)}
    topes = []
    E0 = loops(cocircuits)

    while covectors_new != set():
        Y = covectors_new.pop()
        for X in cocircuits:
            if not X <= Y:  # otherwise Z = X.compose(Y) = Y in F
                Z = X.compose(Y)
                if Z not in covectors:
                    covectors.add(Z)
                    if Z.zero_support() == E0:
                        topes.append(Z)
                    else:
                        covectors_new.add(Z)
    return topes


def lower_faces(covectors):
    r"""
    Compute a list of lower faces.

    INPUT:

    - ``covectors`` -- a list of all covectors with same rank ``r`` of an oriented matroid.

    OUTPUT:

    Returns a list of covectors of rank ``r-1`` of the oriented matroid.

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
        raise ValueError('List is empty.')
    n = covectors[0].length()
    W_ = set()
    for Wj in classes_same_support(covectors):
        PC = parallel_classes(Wj)
        for X in Wj:
            for D in PC:
                for i in D:
                    if X[i] != 0:  # hence X_D != 0
                        if X.reverse_signs_in(D) in Wj:
                            W_.add(sign_vector([0 if i in D else X[i] for i in range(n)]))
                        break
    return list(W_)


def face_enumeration(covectors):
    r"""
    Compute all covectors with less rank than the given list of covectors.

    INPUT:

    - ``covectors`` -- a list of all covectors of same rank ``r`` of an oriented matroid.

    OUTPUT:

    Returns a list of lists. Every list consists of all covectors of the same rank
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

        sage: from sign_vectors.oriented_matroids import topes_from_matrix, face_enumeration
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        [(--+), (+++), (---), (+-+), (-+-), (++-)]
        sage: face_enumeration(tA)
        [[(000)],
         [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
         [(--+), (+++), (---), (+-+), (-+-), (++-)]]
    """
    if not covectors:
        raise ValueError('List is empty.')
    faces = [covectors]

    while faces[0] != [0] and faces[0] != []:
        faces.insert(0, lower_faces(faces[0]))
    return faces


def topes_from_matrix(M, kernel=True):
    r"""
    Return a list of topes of the oriented matroid corresponding to the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``kernel`` is true, returns a list of topes determined by the kernel of
      the matrix ``M``.

    - If ``kernel`` is false, returns a list of topes determined by the row space
      of the matrix ``M``.

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import topes_from_matrix
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: topes_from_matrix(A)
        [(--+), (+++), (---), (+-+), (-+-), (++-)]

    TESTS::

        sage: M = zero_matrix(1, 3)
        sage: topes_from_matrix(M, kernel=False)
        [(000)]
    """
    cc = cocircuits_from_matrix(M, kernel=kernel)
    if cc == []:
        return [zero_sign_vector(M.ncols())]
    return topes_from_cocircuits(cc)


def covectors_from_topes(topes, separate=False):
    r"""
    Compute all covectors from a list of topes.

    INPUT:

    - ``topes`` -- a list of topes.

    - ``separate`` -- a boolean (default: ``False``)

    OUTPUT:

    The list of covectors of the corresponding oriented matroid.

    - If ``separate`` is false, returns a list of covectors. The covectors are
      sorted by rank. (default)

    - If ``separate`` is true, returns a list of lists of covectors, separated
      by their rank.

    .. SEEALSO::

        :func:`~face_enumeration`

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import topes_from_matrix, covectors_from_topes
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        [(--+), (+++), (---), (+-+), (-+-), (++-)]
        sage: covectors_from_topes(tA)
        [(000),
         (-0-),
         (0-+),
         (--0),
         (+0+),
         (0+-),
         (++0),
         (--+),
         (+++),
         (---),
         (+-+),
         (-+-),
         (++-)]
        sage: covectors_from_topes(tA, separate=True)
        [[(000)],
         [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
         [(--+), (+++), (---), (+-+), (-+-), (++-)]]
    """
    if separate:
        return face_enumeration(topes)
    else:
        return flatten(face_enumeration(topes))


def cocircuits_from_topes(topes):
    r"""
    Compute all cocircuits from a list of topes.

    INPUT:

    - ``topes`` -- a list of topes.

    OUTPUT:

    A list of cocircuits of the corresponding oriented matroid.

    EXAMPLES:

    We define some matrix and compute the topes of the corresponding
    oriented matroid::

        sage: from sign_vectors.oriented_matroids import topes_from_matrix, cocircuits_from_topes
        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: tA = topes_from_matrix(A)
        sage: tA
        [(--+), (+++), (---), (+-+), (-+-), (++-)]
        sage: cocircuits_from_topes(tA)
        [(-0-), (0-+), (--0), (+0+), (0+-), (++0)]
    """
    return face_enumeration(topes)[1]


def covectors_from_matrix(M, kernel=True, algorithm=None, separate=False):
    r"""
    Return a list of covectors of the oriented matroid corresponding to the matrix ``A``.

    INPUT:

    - ``M`` -- a matrix.

    - ``kernel`` -- a boolean (default: ``True``)

    - ``algorithm`` -- either ``None`` (default), ``"face_enumeration"`` or ``"fe"``

    - ``separate`` -- a boolean (default: ``False``)

    OUTPUT:

    Returns the list of covectors of an oriented matroid corresponding to the
    matrix ``M``.

    - If ``kernel`` is true, the returned covectors will be determined by the
      kernel of the matrix ``M``.

    - If ``kernel`` is false, the returned covectors will be determined by the
      row space of the matrix ``M``.

    - If ``algorithm`` is ``"face_enumeration"`` or the shortcut ``"fe"``,
      applies the algorithm face enumeration.

    - If ``separate`` is true, returns a list of lists of covectors, separated
      by their rank by applying the algorithm face enumeration.

    - If ``separate`` is false, returns a list of covectors.
      This list is sorted by rank if the algorithm face enumeration is used.

    .. SEEALSO::

        :func:`~face_enumeration`

    EXAMPLES::

        sage: from sign_vectors.oriented_matroids import covectors_from_matrix

    We define some matrix::

        sage: A = matrix([[2, -1, -1]])
        sage: A
        [ 2 -1 -1]
        sage: covectors_from_matrix(A)
        [(---),
         (--+),
         (-0-),
         (000),
         (0-+),
         (--0),
         (+-+),
         (+0+),
         (+++),
         (++-),
         (-+-),
         (0+-),
         (++0)]
        sage: covectors_from_matrix(A, separate=True)
        [[(000)],
         [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
         [(--+), (+++), (---), (+-+), (-+-), (++-)]]
        sage: covectors_from_matrix(A, algorithm="face_enumeration", separate=True)
        [[(000)],
         [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
         [(--+), (+++), (---), (+-+), (-+-), (++-)]]
        sage: covectors_from_matrix(A, algorithm="fe", separate=True)
        [[(000)],
         [(-0-), (0-+), (--0), (+0+), (0+-), (++0)],
         [(--+), (+++), (---), (+-+), (-+-), (++-)]]

    TESTS::

        sage: covectors_from_matrix(zero_matrix(1, 4), kernel=False, separate=False)
        [(0000)]
        sage: covectors_from_matrix(zero_matrix(1, 4), kernel=False, separate=True)
        [[(0000)]]
    """
    if algorithm is None:
        if separate:
            algorithm = 'face_enumeration'
        else:
            cc = cocircuits_from_matrix(M, kernel=kernel)
            if cc == []:
                return [zero_sign_vector(M.ncols())]
            return covectors_from_cocircuits(cc)
    if algorithm in ['face_enumeration', 'fe']:
        return covectors_from_topes(topes_from_matrix(M, kernel=kernel), separate=separate)
    else:
        raise ValueError("no algorithm '" + algorithm + "'")

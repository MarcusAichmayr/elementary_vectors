r"""Utility functions and other useful functions for working with oriented matroids."""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.modules.free_module_element import vector
from sign_vectors import SignVector, sign_vector


def _subvector(F, R):
    r"""
    Return a function that returns a sign vector or vector consisting of entries not in ``R``. Used by ``contraction`` and ``deletion``.

    INPUT:

    - ``F`` -- a list of vectors or sign vectors

    - ``R`` -- a list of indices

    EXAMPLES::

        sage: from sign_vectors.utility import _subvector
        sage: from sign_vectors import sign_vector
        sage: W = [sign_vector("++0"), sign_vector("-00"), sign_vector("00+")]
        sage: W
        [(++0), (-00), (00+)]
        sage: f = _subvector(W, [1])
        sage: f(sign_vector("-+0"))
        (-0)
        sage: l = [vector([0,0,1]), vector([0,2,1]), vector([-1,0,1])]
        sage: f = _subvector(l, [1])
        sage: f(vector([1,2,3]))
        (1, 3)
    """
    if F == []:
        raise ValueError('List is empty.')
    n = len(list(F[0]))
    S = [e for e in range(n) if e not in R]  # S = E\R

    if isinstance(F[0], SignVector):
        def vec(v):
            return sign_vector(v.list_from_positions(S))
    else:
        def vec(v):
            return vector(v.list_from_positions(S))
    return vec


def loops(W):
    r"""
    Compute the list of loops of a given list of sign vectors (or real vectors) ``W``.

    .. NOTE::

        A loop is a component where every sign vector of ``W`` is zero.

    EXAMPLES::

        sage: from sign_vectors.utility import loops
        sage: from sign_vectors import sign_vector
        sage: loops([sign_vector([0,1,0]), sign_vector([-1,0,0])])
        [2]
        sage: loops([sign_vector([1,0,0]), sign_vector([-1,0,0])])
        [1, 2]

    Also works for real vectors::

        sage: loops([vector([5,0,0,0]), vector([2,0,-3,0])])
        [1, 3]
    """
    if not W:
        raise ValueError('List is empty.')
    n = W[0].length()

    L = []
    for e in range(n):
        val = True
        for X in W:
            if X[e] != 0:
                val = False
                break
        if val is True:
            L.append(e)

    return L


def is_parallel(W, e, f, return_ratio=False):
    r"""
    Determine whether two elements ``e, f`` are parallel for each vector of ``W``. This also works for a set of sign vectors.

    INPUT:

    - ``W`` -- a list of vectors or sign vectors of length ``n``

    - ``e`` -- an integer with :math:`0 \leq e \leq n-1`

    - ``f`` -- an integer with :math:`0 \leq e \leq n-1`

    - ``return_ratio`` -- a boolean (default: False)

    OUTPUT:

    Returns a boolean.
    If ``return_ratio`` is true, a list consisting of the boolean and the ratio will be returned instead.

    .. NOTE::

        The elements ``e`` and ``f`` are parallel if there exists a ratio ``d`` such that
        ``v[e] = d v[f]`` for each ``v`` in ``W``.

    EXAMPLES::

        sage: from sign_vectors.utility import is_parallel
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: is_parallel(L, 0, 1)
        False
        sage: is_parallel(L, 1, 2)
        False
        sage: is_parallel(L, 1, 3)
        True

    Now, we consider some real vectors::

        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: is_parallel(L, 0, 2)
        True
        sage: is_parallel(L, 0, 1)
        False
        sage: is_parallel(L, 1, 3)
        False
        sage: is_parallel(L, 4, 5)
        True

    We can also return the ratio of the two components::

        sage: is_parallel(L, 0, 2, return_ratio=True)
        [True, 1/2]
        sage: is_parallel(L, 2, 0, return_ratio=True)
        [True, 2]
        sage: is_parallel(L, 0, 1, return_ratio=True)
        [False, 0]

    Also works for matrices::

        sage: M = matrix([[0,0,1,-1,0],[1,0,0,0,1],[1,1,1,1,1]])
        sage: is_parallel(M,0,4)
        True
        sage: is_parallel(M,0,1)
        False
    """
    d = 0  # will later be set to the ratio of X[e] and X[f]

    if return_ratio:
        ret = [False, 0]
    else:
        ret = False

    for X in W:
        if d == 0:
            if X[f] == 0:
                if X[e] != 0:
                    return ret
            elif X[e] == 0:
                return ret
            else:  # determine ratio
                d = X[e]/X[f]
        else:
            if X[e] != d * X[f]:
                return ret
    if return_ratio:
        return [True, d]
    else:
        return True


def parallel_classes(W, positive_only=False):
    r"""
    Compute the parallel classes of a given set of vectors ``W``. This also works for a set of sign vectors.

    INPUT:

    - ``W`` -- a list of vectors or sign vectors of length ``n``

    - ``positive_only`` -- a boolean (default: False)

    OUTPUT:

    Returns a partition of ``[0, ..., n-1]`` into parallel classes.

    If ``positive_only`` is true, returns a partition of ``[0, ..., n-1]`` into positive parallel classes,
    that is, the ratios of the corresponding classes are non-negative.

    .. NOTE::

        The elements ``e`` and ``f`` are parallel if there exists a ratio ``d`` such that
        ``v[e] = d v[f]`` for each ``v`` in ``W``.

    EXAMPLES::

        sage: from sign_vectors.utility import parallel_classes
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: parallel_classes(L)
        [[0], [1, 3], [2]]
        sage: parallel_classes(L, positive_only=True)
        [[0], [1], [2], [3]]

    Now, we compute the parallel classes of a list of real vectors::

        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: parallel_classes(L)
        [[0, 2], [1], [3], [4, 5]]

    Let us compute the parallel classes of the rows of a matrix::

        sage: M = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: parallel_classes(M)
        [[0, 4], [1], [2, 3]]
        sage: parallel_classes(M, positive_only=True)
        [[0, 4], [1], [2], [3]]
    """
    if not W:
        raise ValueError('List is empty.')
    L = []
    k = W[0].length()
    toCheck = list(range(k))

    if positive_only:
        def is_par(W, e, f):
            val = is_parallel(W, e, f, return_ratio=True)
            return val[1] > 0 if val[0] else False
    else:
        def is_par(W, e, f):
            return is_parallel(W, e, f)

    while len(toCheck) > 0:
        e = toCheck.pop(0)
        l = [e]
        for f in toCheck[:]:  # find parallel class ``l`` of ``e``
            if is_par(W, e, f):
                l.append(f)
                toCheck.remove(f)
        L.append(l)
    return L


def positive_parallel_classes(W):
    r"""
    Compute the positive parallel classes of a given set of vectors ``W``. This also works for a set of sign vectors.

    .. SEEALSO::

        :func:`~parallel_classes`

    EXAMPLES::

        sage: from sign_vectors.utility import positive_parallel_classes
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("--0+"), sign_vector("00+0")]
        sage: L
        [(++0-), (--0+), (00+0)]
        sage: positive_parallel_classes(L)
        [[0, 1], [2], [3]]

    Now, we compute the positive parallel classes of a list of real vectors::

        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: positive_parallel_classes(L)
        [[0, 2], [1], [3], [4], [5]]

    Let us compute the positive parallel classes of the rows of a matrix::

        sage: M = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: positive_parallel_classes(M)
        [[0, 4], [1], [2], [3]]
    """
    return parallel_classes(W, positive_only=True)


def classes_same_support(W):
    r"""
    Compute the classes with same support of a given list of sign vectors. Also works for a list of real vectors.

    INPUT:

    - ``W`` -- a list of vectors of length ``n``.

    EXAMPLES::

        sage: from sign_vectors.utility import classes_same_support
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: classes_same_support(L)
        [[(++0-), (+-0+)], [(-0+0)]]
        sage: classes_same_support([vector([1,1,0,0]), vector([2,-3,0,0]), vector([0,1,0,0])])
        [[(1, 1, 0, 0), (2, -3, 0, 0)], [(0, 1, 0, 0)]]
    """
    L = dict()
    for X in W:
        s = tuple(X.support())  # tuples are hashable
        if s not in L.keys():
            L[s] = [X]
        else:
            L[s].append(X)
    return list(L.values())


def adjacent(X, Y, S):
    r"""
    Return whether the sign vectors ``X`` and ``Y`` are adjacent over the set of sign vectors ``S``.

    INPUT:

    - ``X`` -- a sign vector

    - ``Y`` -- a sign vector

    - ``S`` -- a list of sign vectors

    OUTPUT:
    a boolean

    .. NOTE::

        define adjacent here TODO

    EXAMPLES:

    We consider the following matrix::

        sage: M = matrix([[1,2,0],[0,1,-1]])
        sage: M
        [ 1  2  0]
        [ 0  1 -1]

    By using the function :func:`sign_vectors.oriented_matroids.cocircuits_from_matrix`, we can compute the corresponding cocircuits::

        sage: from sign_vectors.oriented_matroids import *
        sage: cc = cocircuits_from_matrix(M)
        sage: cc
        [(--0), (++0), (-0-), (+0+), (0-+), (0+-)]

    The two sign vectors ``X = (++0)`` and ``Y = (+0+)`` are harmonious::

        sage: X = sign_vector('++0')
        sage: X
        (++0)
        sage: Y = sign_vector('+0+')
        sage: Y
        (+0+)
        sage: X.is_harmonious(Y)
        True

    Furthermore, the only cocircuits lying under the composition of :math:`X` and :math:`Y`,
    that is, cocircuits :math:`Z` satisfying :math:`Z < (+++) = X \circ Y`,
    are :math:`X` and :math:`Y`.
    Hence, those two sign vectors are adjacent::

        sage: from sign_vectors.utility import adjacent
        sage: adjacent(X, Y, cc)
        True

    Conversely, :math:`Y = (+0+)` and :math:`Z = (0+-)` are not adjacent since
    :math:`(++0) < (++-) = Y \circ Z`::

        sage: Z = sign_vector('0+-')
        sage: Z
        (0+-)
        sage: adjacent(Y, Z, cc)
        False
    """
    XY = X & Y
    return not any(Z < XY for Z in S if Z != X and Z != Y)

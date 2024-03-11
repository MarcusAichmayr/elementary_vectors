"""Utility functions and other useful functions for working with oriented matroids"""

#############################################################################
#  Copyright (C) 2024                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.modules.free_module_element import vector
from sign_vectors import SignVector, sign_vector


def loops(iterable):
    r"""
    Compute the loops of sign vectors or vectors.

    .. NOTE::

        A loop is a component where every element is zero.

    EXAMPLES::

        sage: from sign_vectors.utility import loops
        sage: from sign_vectors import sign_vector
        sage: loops([sign_vector([0, 1, 0]), sign_vector([-1, 0, 0])])
        [2]
        sage: loops([sign_vector([1, 0, 0]), sign_vector([-1, 0, 0])])
        [1, 2]

    Also works for real vectors::

        sage: loops([vector([5, 0, 0, 0]), vector([2, 0, -3, 0])])
        [1, 3]
    """
    if not iterable:
        raise ValueError('Iterable is empty.')
    for _ in iterable:
        length = _.length()
        break

    return [e for e in range(length) if all(element[e] == 0 for element in iterable)]


def is_parallel(iterable, component1, component2, return_ratio=False):
    r"""
    Determine whether two components of sign vectors or vectors are parallel.

    INPUT:

    - ``iterable`` -- a list of sign vectors or vectors of length ``n``

    - ``component1`` -- an integer with ``0 <= component1 < n``

    - ``component2`` -- an integer with ``0  <= component2 < n``

    - ``return_ratio`` -- a boolean

    OUTPUT:

    Returns a boolean.
    If ``return_ratio`` is true, a boolean and the ratio will be returned instead.

    .. NOTE::

        The elements ``component1`` and ``component2`` are parallel if there exists a ratio ``d`` such that
        ``v[component1] = d v[component2]`` for each ``v`` in ``iterable``.

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

        sage: L = [vector([1, 1, 2, 3, 0, 0]), vector([-2, 1, -4, 3, 3, -17]), vector([0, 1, 0, 1, 0, 0])]
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
        (True, 1/2)
        sage: is_parallel(L, 2, 0, return_ratio=True)
        (True, 2)
        sage: is_parallel(L, 0, 1, return_ratio=True)
        (False, 0)

    Also works for matrices::

        sage: M = matrix([[0, 0, 1, -1, 0], [1, 0, 0, 0, 1], [1, 1, 1, 1, 1]])
        sage: is_parallel(M, 0, 4)
        True
        sage: is_parallel(M, 0, 1)
        False
    """
    ratio = 0

    if return_ratio:
        false_return = False, 0
    else:
        false_return = False

    for element in iterable:
        if ratio == 0:
            if element[component2] == 0:
                if element[component1] != 0:
                    return false_return
            elif element[component1] == 0:
                return false_return
            else:
                ratio = element[component1] / element[component2]
        else:
            if element[component1] != ratio * element[component2]:
                return false_return
    if return_ratio:
        return True, ratio
    return True


def parallel_classes(iterable, positive_only=False):
    r"""
    Compute the parallel classes of given sign vectors or vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors with same length

    - ``positive_only`` -- a boolean (default: False)

    OUTPUT:

    Returns a partition of ``[0, ..., n - 1]`` into parallel classes.

    If ``positive_only`` is true, returns a partition of ``[0, ..., n - 1]`` into positive parallel classes,
    that is, the ratios of the corresponding classes are non-negative.

    .. NOTE::

        The elements ``component1`` and ``component2`` are parallel if there exists a ratio ``d`` such that
        ``v[component1] = d v[component2]`` for each ``v`` in ``iterable``.

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

        sage: L = [vector([1, 1, 2, 3, 0, 0]), vector([-2, 1, -4, 3, 3, -17]), vector([0, 1, 0, 1, 0, 0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: parallel_classes(L)
        [[0, 2], [1], [3], [4, 5]]

    Let us compute the parallel classes of the rows of a matrix::

        sage: M = matrix([[0, 0, 1, -2, 0], [1, 0, 0, 0, 1], [1, 1, -3, 6, 1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: parallel_classes(M)
        [[0, 4], [1], [2, 3]]
        sage: parallel_classes(M, positive_only=True)
        [[0, 4], [1], [2], [3]]
    """
    if not iterable:
        raise ValueError('Iterable is empty.')
    output = []
    indices_to_check = list(range(iterable[0].length()))

    if positive_only:
        def is_par(iterable, component1, component2):
            value = is_parallel(iterable, component1, component2, return_ratio=True)
            return value[1] > 0 if value[0] else False
    else:
        def is_par(iterable, component1, component2):
            return is_parallel(iterable, component1, component2)

    while len(indices_to_check) > 0:
        component1 = indices_to_check.pop(0)
        parallel_class = [component1]
        for component2 in indices_to_check[:]:
            if is_par(iterable, component1, component2):
                parallel_class.append(component2)
                indices_to_check.remove(component2)
        output.append(parallel_class)
    return output


def positive_parallel_classes(iterable):
    r"""
    Compute the positive parallel classes of given sign vectors or vectors.

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

        sage: L = [vector([1, 1, 2, 3, 0, 0]), vector([-2, 1, -4, 3, 3, -17]), vector([0, 1, 0, 1, 0, 0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: positive_parallel_classes(L)
        [[0, 2], [1], [3], [4], [5]]

    Let us compute the positive parallel classes of the rows of a matrix::

        sage: M = matrix([[0, 0, 1, -2, 0], [1, 0, 0, 0, 1], [1, 1, -3, 6, 1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: positive_parallel_classes(M)
        [[0, 4], [1], [2], [3]]
    """
    return parallel_classes(iterable, positive_only=True)


def classes_same_support(iterable):
    r"""
    Compute the classes with same support of given sign vectors or vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors with same length

    EXAMPLES::

        sage: from sign_vectors.utility import classes_same_support
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: classes_same_support(L)
        [[(++0-), (+-0+)], [(-0+0)]]
        sage: classes_same_support([vector([1, 1, 0, 0]), vector([2, -3, 0, 0]), vector([0, 1, 0, 0])])
        [[(1, 1, 0, 0), (2, -3, 0, 0)], [(0, 1, 0, 0)]]
    """
    output = {}
    for element in iterable:
        support = tuple(element.support())  # tuples are hashable
        if support not in output.keys():
            output[support] = [element]
        else:
            output[support].append(element)
    return list(output.values())


def adjacent(element1, element2, iterable):
    r"""
    Return whether two sign vectors are adjacent over given sign vectors.

    INPUT:

    - ``element1`` -- a sign vector

    - ``element2`` -- a sign vector

    - ``iterable`` -- an iterable of sign vectors

    OUTPUT:
    a boolean

    .. NOTE::

        define adjacent here TODO

    EXAMPLES:

    We consider the following matrix::

        sage: M = matrix([[1, 2, 0], [0, 1, -1]])
        sage: M
        [ 1  2  0]
        [ 0  1 -1]

    By using the function :func:`sign_vectors.oriented_matroids.cocircuits_from_matrix`, we can compute the corresponding cocircuits::

        sage: from sign_vectors.oriented_matroids import *
        sage: cc = cocircuits_from_matrix(M, kernel=False)
        sage: cc
        {(0-+), (+0+), (--0), (-0-), (0+-), (++0)}

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
    composition = element1 & element2
    return not any(Z < composition for Z in iterable if not Z in [element1, element2])


def exclude_indices(vectors, indices):
    r"""
    Return a function that returns a sign vector or vector with entries not in given indices.
    
    INPUT:

    - ``vectors`` -- a list of sign vectors or vectors

    - ``indices`` -- a list of indices

    EXAMPLES::

        sage: from sign_vectors.utility import exclude_indices
        sage: from sign_vectors import sign_vector
        sage: W = [sign_vector("++0"), sign_vector("-00"), sign_vector("00+")]
        sage: W
        [(++0), (-00), (00+)]
        sage: f = exclude_indices(W, [1])
        sage: f(sign_vector("-+0"))
        (-0)
        sage: l = [vector([0, 0, 1]), vector([0, 2, 1]), vector([-1, 0, 1])]
        sage: f = exclude_indices(l, [1])
        sage: f(vector([1, 2, 3]))
        (1, 3)
    """
    if not vectors:
        raise ValueError('List is empty.')
    length = len(list(vectors[0]))
    other_indices = [e for e in range(length) if not e in indices]

    if isinstance(vectors[0], SignVector):
        def vec(iterable):
            return sign_vector(iterable.list_from_positions(other_indices))
    else:
        def vec(iterable):
            iterable = vector(iterable.list_from_positions(other_indices))
            iterable.set_immutable()
            return iterable
    return vec

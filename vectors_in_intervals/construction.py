"""
Constructing vectors with components in intervals

The core function of this module is :func:`~construct_vector`.
Most other functions here are auxiliary functions but can be used for special cases.
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

from sage.misc.mrange import cartesian_product_iterator
from sage.modules.free_module_element import vector, zero_vector

from sign_vectors import sign_vector, SignVector
from elementary_vectors.functions import ElementaryVectors
from . import Intervals


def vector_from_sign_vector(data, sv: SignVector) -> vector:
    r"""
    Find a vector in the row space of a matrix that has given signs.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                elementary vectors of length ``n``

    - ``sv`` -- a sign vector of length ``n``

    OUTPUT:
    Return a conformal sum of elementary vectors that lies in the given subspace.

    If ``data`` is a matrix, the elementary vectors in the kernel of this matrix are used for the result.
    If ``data`` is a list of elementary vectors, those are used.

    .. NOTE::

        A ``ValueError`` is raised if no solution exists.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: from elementary_vectors import *
        sage: from sign_vectors import *
        sage: M = matrix([[1, 0, 2, 0], [0, 1, 1, 0], [0, 0, 0, 1]])
        sage: vector_from_sign_vector(M, zero_sign_vector(4))
        (0, 0, 0, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-+0"))
        (2, -2, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+0+0"))
        (1, 0, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: evs = elementary_vectors(M, dual=False)
        sage: vector_from_sign_vector(evs, sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: vector_from_sign_vector(M, sign_vector("+0-0"))
        Traceback (most recent call last):
        ...
        ValueError: Cannot find vector corresponding to given sign vector.
        sage: vector_from_sign_vector([], zero_sign_vector(4))
        (0, 0, 0, 0)
    """
    try:
        return vector_between_sign_vectors(data, sv, sv)
    except ValueError as exc:
        raise ValueError("Cannot find vector corresponding to given sign vector.") from exc


def vector_between_sign_vectors(data, lower: SignVector, upper: SignVector) -> vector:
    r"""
    Find a vector in the row space of a matrix that has given signs.

    The resulting vector ``v`` satisfies ``lower <= sign(v) <= upper``.

    .. SEEALSO::

        :func:`~vector_from_sign_vector`

    TESTS::

        sage: from vectors_in_intervals import *
        sage: from sign_vectors import *
        sage: M = matrix([[-3, -8, 0, 0, 1, -1], [10, -1, 2, 1, -7, 0], [1, 0, 0, -1, -3, 3]])
        sage: lower = sign_vector('000+00')
        sage: lower
        (000+00)
        sage: upper = sign_vector('++0+0+')
        sage: upper
        (++0+0+)

    We demonstrate that we cannot just use evs with indices in supp X::

        sage: M = matrix([[1, 1, 1, 0, 0], [0, 0, 0, 1, 1]])
        sage: X = sign_vector("+++00")
        sage: vector_between_sign_vectors(M, X, X)
        (1, 1, 1, 0, 0)
    """
    if isinstance(data, list):
        evs = data
        try:
            result = data[0].parent().zero_vector()
        except IndexError:
            result = zero_vector(lower.length())
    elif isinstance(data, Generator):
        evs = data
        result = zero_vector(lower.length())
    else:
        evs_object = ElementaryVectors(data)
        # evs_object.set_combinations_dual(Combinations(upper.support(), evs_object.length - evs_object.rank + 1))
        evs = evs_object.generator(dual=False)
        result = zero_vector(data.base_ring(), lower.length())

    if sign_vector(result) >= lower:
        return result
    for v in evs:
        for w in [v, -v]:
            if sign_vector(w) <= upper:
                result += w
                if sign_vector(result) >= lower:
                    return result
                break

    raise ValueError("Cannot find vector corresponding to given sign vectors.")


def sign_vectors_in_intervals(intervals: Intervals, generator: bool = False) -> list[SignVector] | Generator[SignVector]:
    r"""
    Compute all sign vectors that correspond to a vector with components in given intervals.

    INPUT:

    - ``intervals`` -- a list of intervals

    - ``generator`` -- a boolean (default: ``False``)

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: intervals = Intervals.from_bounds([-1, 1], [0, 1])
        sage: sign_vectors_in_intervals(intervals)
        [(0+), (-+)]
        sage: intervals = Intervals.from_bounds([-1, -2], [0, 1])
        sage: sign_vectors_in_intervals(intervals)
        [(00), (0+), (0-), (-0), (-+), (--)]
        sage: intervals = Intervals.from_bounds([-1, -1, 0], [0, 5, 0])
        sage: sign_vectors_in_intervals(intervals)
        [(000), (0+0), (0-0), (-00), (-+0), (--0)]
        sage: intervals = Intervals.from_bounds([-1, -1, -1], [0, 1, 0], False, False)
        sage: sign_vectors_in_intervals(intervals)
        [(-0-), (-+-), (---)]
        sage: intervals = Intervals.from_bounds([-1, 0], [1, 0], False, False)
        sage: sign_vectors_in_intervals(intervals)
        []
        sage: intervals = Intervals.from_bounds([], [])
        sage: sign_vectors_in_intervals(intervals)
        []
    """
    list_of_signs = []
    if intervals.is_empty():
        if generator:
            def empty():
                yield from ()
            return empty()
        return []
    for interval in intervals:
        available_signs = []
        if 0 in interval:
            available_signs.append(0)
        if interval.supremum() > 0:
            available_signs.append(1)
        if interval.infimum() < 0:
            available_signs.append(-1)
        list_of_signs.append(available_signs)

    if generator:
        return (
            sign_vector(signs) for signs in cartesian_product_iterator(list_of_signs)
        )
    return [sign_vector(signs) for signs in cartesian_product_iterator(list_of_signs)]

"""Existence of vectors with components in intervals"""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.sets.real_set import RealSet

from elementary_vectors import elementary_vectors


def exists_orthogonal_vector(v, intervals: list[RealSet]) -> bool:
    r"""
    Check whether an orthogonal vector exists with components in given intervals.

    INPUT:

    - ``v`` -- a vector of length ``n``

    - ``intervals`` -- a list of ``n`` intervals

    OUTPUT:

    Return whether there exists an orthogonal vector to ``v``
    such that the components of this vector lie in ``intervals``.

    EXAMPLES:

    We define several lists of intervals and vectors
    and apply the function::

        sage: from vectors_in_intervals import *
        sage: I = intervals_from_bounds([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1, 1, 1])
        sage: exists_orthogonal_vector(v, I)
        True
        sage: v = vector([1, 1, -1])
        sage: exists_orthogonal_vector(v, I)
        False

    Next, we consider open intervals::

        sage: I = intervals_from_bounds([0, 1, -1], [1, 2, 1], False, [True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1, 0, 1])
        sage: exists_orthogonal_vector(v, I)
        True

    We can also work with unbounded intervals::

        sage: I = intervals_from_bounds([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1, 1, -1])
        sage: exists_orthogonal_vector(v, I)
        True

    TESTS::

        sage: I = intervals_from_bounds([-oo, 0], [oo, 0], False, False)
        sage: I
        [(-oo, +oo), {}]
        sage: v = vector([1, 0])
        sage: exists_orthogonal_vector(v, I)
        False
        sage: v = vector([1])
        sage: exists_orthogonal_vector(v, I)
        Traceback (most recent call last):
        ...
        ValueError: Lengths of vector and intervals do not match.
    """
    if len(v) != len(intervals):
        raise ValueError("Lengths of vector and intervals do not match.")

    lower_product = 0
    upper_product = 0
    lower_product_attainable = True
    upper_product_attainable = True

    for entry, interval in zip(v, intervals):
        if interval.is_empty():
            return False
        if not entry:
            continue
        bound1 = interval.inf() if entry > 0 else interval.sup()
        bound2 = interval.sup() if entry > 0 else interval.inf()
        lower_product += entry * bound1
        upper_product += entry * bound2
        lower_product_attainable &= bound1 in interval
        upper_product_attainable &= bound2 in interval

    if lower_product > 0:
        return False
    if upper_product < 0:
        return False
    if lower_product == 0 and not lower_product_attainable:
        return False
    if upper_product == 0 and not upper_product_attainable:
        return False
    return True


def exists_vector(data, intervals: list[RealSet], certify: bool = False) -> bool:
    r"""
    Check whether a vector exists in a given vector space with components in given intervals.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                  elementary vectors of length ``n``

    - ``intervals`` -- a list of ``n`` intervals

    - ``certify`` -- a boolean (default: ``False``)

    OUTPUT:

    Return whether there exists a vector in a given vector space
    such that the components lie in specified intervals using elementary vectors.

    - If ``data`` is a matrix, check if a vector in the row space lies in the intervals.

    - If ``data`` is a list of elementary vectors, check if a vector exists orthogonal to those.

    - If ``certify`` is true and no vector exists, an elementary vector is returned as certificate.

    ALGORITHM:

    The underlying algorithm is based on Minty's Lemma. (see [Min74]_)

    .. [Min74] Minty, G. J.:
       "A 'from scratch' proof of a theorem of Rockafellar and Fulkerson".
       In: Mathematical Programming 7 (1974), pp. 368-375.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: M = matrix([1, 1, 0])
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    First, we consider closed intervals::

        sage: I = intervals_from_bounds(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: exists_vector(M, I)
        True

    Next, we take open intervals::

        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, False, False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: exists_vector(M, I)
        False

    Since no vector exists, there is an elementary vector certifying this.
    To find one, we pass ``certify=True``::

        sage: exists_vector(M, I, certify=True)
        (1, -1, 0)

    Mixed intervals are also possible::

        sage: lower_bounds_closed = [True, True, False]
        sage: upper_bounds_closed = [False, True, True]
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: I
        [[2, 5), [5, 6], (-1, 1]]
        sage: exists_vector(M, I)
        False

    Finally, we consider unbounded intervals::

        sage: M = matrix([[1, 0, 1, 0], [0, 1, 1, 1]])
        sage: lower_bounds = [2, 5, 0, -oo]
        sage: upper_bounds = [5, oo, 8, 5]
        sage: lower_bounds_closed = [True, True, False, False]
        sage: upper_bounds_closed = [False, False, False, True]
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: exists_vector(M, I)
        True

    TESTS::

        sage: I = intervals_from_bounds([0, 0], [1, 0], False)
        sage: I
        [(0, 1], {}]
        sage: exists_vector([], I)
        False
        sage: M = random_matrix(QQ, 0, 2)
        sage: exists_vector(M, I)
        False
        sage: M = random_matrix(QQ, 0, 1)
        sage: exists_vector(M, I)
        Traceback (most recent call last):
        ...
        ValueError: Number of columns of matrix ``data`` and length of ``intervals`` do not match.
    """
    if hasattr(data, "ncols") and data.ncols() != len(intervals):
        raise ValueError(
            "Number of columns of matrix ``data`` and length of ``intervals`` do not match."
        )
    if any(interval.is_empty() for interval in intervals):
        return False
    evs = data if isinstance(data, list) else elementary_vectors(data, generator=True)
    for v in evs:
        if not exists_orthogonal_vector(v, intervals):
            return v if certify else False
    return True

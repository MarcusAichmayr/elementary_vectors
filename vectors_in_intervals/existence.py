"""Existence of vectors with components in intervals."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from elementary_vectors import elementary_vectors


def exists_orthogonal_vector(v, intervals):
    r"""
    Check whether a normal vector exists that lies in the specified intervals.

    INPUT:

    - ``v`` -- a vector of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    OUTPUT:

    Return whether there exists a vector ``z`` such that the scalar product of ``z`` and ``v`` is zero
    and each component of ``z`` lies in the respective interval of the list ``intervals``.
    Raises a ``ValueError`` if the lengths of ``v`` and ``intervals`` are different.

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

    We can even consider unbounded intervals::

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
        ValueError: Lengths of ``v`` and ``intervals`` are different!
    """
    if len(v) != len(intervals):
        raise ValueError("Lengths of ``v`` and ``intervals`` are different!")

    lower_bound = 0
    upper_bound = 0
    lower_bound_zero = True
    upper_bound_zero = True

    for entry, interval in zip(v, intervals):
        if interval.is_empty():
            return False
        if entry != 0:
            if interval.is_finite():  # interval consists of one point
                lower_bound += entry * interval.an_element()
                upper_bound += entry * interval.an_element()
            elif entry > 0:
                lower_bound += entry * interval.inf()
                upper_bound += entry * interval.sup()
                lower_bound_zero &= interval.inf() in interval
                upper_bound_zero &= interval.sup() in interval
            else:  # entry < 0
                lower_bound += entry * interval.sup()
                upper_bound += entry * interval.inf()
                lower_bound_zero &= interval.sup() in interval
                upper_bound_zero &= interval.inf() in interval

    if lower_bound > 0:
        return False
    if upper_bound < 0:
        return False
    if lower_bound == 0 and not lower_bound_zero:
        return False
    if upper_bound == 0 and not upper_bound_zero:
        return False
    return True


def exists_vector(data, intervals, certify=False):
    r"""
    Return whether a vector exists in a given vector space such that the components lie in the specified intervals.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                  elementary vectors of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)
    
    - ``certify`` -- a boolean (default: ``False``)

    OUTPUT:

    Return whether there exists a vector in a vector space
    such that the components lie in specified intervals using elementary vectors.

    - If ``data`` is a matrix, check if a vector in the row space of this matrix lies in the intervals.

    - If ``data`` is a list of elementary vectors, check if a vector exists orthogonal to those elementary vectors.

    - If ``certify`` is true and no vector exists, an elementary vector certifying non-existence is returned.

    ALGORITHM:

    The underlying algorithm is based on Minty's Lemma. (see [Min74]_)

    .. [Min74] Minty, G. J.:
       „A `from scratch` proof of a theorem of Rockafellar and Fulkerson“.
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
        ValueError: Number of columns of matrix ``data`` and length of ``intervals`` are different!
    """
    if hasattr(data, "ncols"):
        if data.ncols() != len(intervals):
            raise ValueError("Number of columns of matrix ``data`` and length of ``intervals`` are different!")
        for interval in intervals:
            if interval.is_empty():
                return False
        evs = elementary_vectors(data, generator=True)
    else:
        for interval in intervals:
            if interval.is_empty():
                return False
        evs = data

    for v in evs:
        if not exists_orthogonal_vector(v, intervals):
            return v if certify else False
    return True

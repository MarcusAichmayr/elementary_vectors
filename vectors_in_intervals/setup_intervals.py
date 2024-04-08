"""Setting up lists of intervals"""

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
from sage.rings.infinity import Infinity
from sage.rings.rational_field import QQ

from .utility import interval_from_bounds, random_interval


def intervals_from_bounds(
    lower_bounds: list,
    upper_bounds: list,
    lower_bounds_closed: bool = True,
    upper_bounds_closed: bool = True,
) -> list[RealSet]:
    r"""
    Construct a list of intervals from lists of bounds.

    INPUT:

    - ``lower_bounds`` -- a list of real values and infinity of length ``n``

    - ``upper_bounds`` -- a list of real values and infinity of length ``n``

    - ``lower_bounds_closed`` -- a boolean (default: ``True``) or a list of booleans of length ``n``

    - ``upper_bounds_closed`` -- a boolean (default: ``True``) or a list of booleans of length ``n``

    OUTPUT:

    A list of ``RealSet`` objects of length ``n``.

    - ``lower_bounds`` and ``upper_bounds`` are the lower and upper interval bounds.
      If ``lower_bounds[i] > upper_bounds[i]``, those elements will be exchanged.

    - ``lower_bounds_closed`` and ``upper_bounds_closed`` determine the intervals.

    - The lower (or upper) interval bounds of the ``i``-th interval is

        - closed if ``lower_bounds_closed[i]`` (or ``upper_bounds_closed[i]``) is ``True`` (default).

        - open if ``lower_bounds_closed[i]`` (or ``upper_bounds_closed[i]``) is ``False``.

    - If ``lower_bounds_closed`` (or ``upper_bounds_closed``) is a boolean,
      then all lower or upper interval bounds
      are considered closed if ``True`` (default) and open if ``False``.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    By default, the intervals are closed::

        sage: intervals_from_bounds(lower_bounds, upper_bounds)
        [[2, 5], [5, 6], [-1, 1]]

    We obtain open intervals if both ``lower_bounds_closed`` and ``upper_bounds_closed`` are false::

        sage: intervals_from_bounds(lower_bounds, upper_bounds, False, False)
        [(2, 5), (5, 6), (-1, 1)]

    Mixed intervals are also possible::

        sage: intervals_from_bounds(lower_bounds, upper_bounds, False)
        [(2, 5], (5, 6], (-1, 1]]
        sage: intervals_from_bounds(lower_bounds, upper_bounds, [False, True, True], [True, True, False])
        [(2, 5], [5, 6], [-1, 1)]

    We can also specify unbounded intervals.
    Note that bounds at infinity are always open::

        sage: lower_bounds = [-oo, 2, -oo]
        sage: upper_bounds = [-5, oo, oo]
        sage: intervals_from_bounds(lower_bounds, upper_bounds)
        [(-oo, -5], [2, +oo), (-oo, +oo)]

    Finite and empty intervals are represented as usual::

        sage: lower_bounds = [0, 2]
        sage: upper_bounds = [0, 2]
        sage: intervals_from_bounds(lower_bounds, upper_bounds, [True, False], True)
        [{0}, {}]

    If the lower bound is greater than the upper bound, those bounds will be exchanged::

        sage: lower_bounds = [1, 4]
        sage: upper_bounds = [2, 3]
        sage: intervals_from_bounds(lower_bounds, upper_bounds, False)
        [(1, 2], (3, 4]]
    """
    length = len(lower_bounds)
    if len(upper_bounds) != length:
        raise ValueError(
            "``upper_bounds`` should be a list of length " + str(length) + "."
        )

    if lower_bounds_closed is True:
        lower_bounds_closed = [True] * length
    elif lower_bounds_closed is False:
        lower_bounds_closed = [False] * length
    elif len(lower_bounds_closed) != length:
        raise ValueError(
            "``lower_bounds_closed`` should be a list of length " + str(length) + "."
        )

    if upper_bounds_closed is True:
        upper_bounds_closed = [True] * length
    elif upper_bounds_closed is False:
        upper_bounds_closed = [False] * length
    elif len(upper_bounds_closed) != length:
        raise ValueError(
            "``upper_bounds_closed`` should be a list of length " + str(length) + "."
        )

    return [
        interval_from_bounds(*bounds)
        for bounds in zip(
            lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed
        )
    ]


def is_vector_in_intervals(v, intervals: list[RealSet]) -> bool:
    r"""
    Check if a vector lies in a list of intervals.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: v = vector([5, 0, -10])
        sage: intervals = intervals_from_bounds([0, 0, -oo], [oo, 0, 0])
        sage: is_vector_in_intervals(v, intervals)
        True
        sage: is_vector_in_intervals(-v, intervals)
        False
    """
    return all(entry in interval for entry, interval in zip(v, intervals))


def random_intervals(
    length: int,
    ring=QQ,
    allow_infinity: bool = True,
    allow_open: bool = True,
    allow_empty: bool = False,
) -> list[RealSet]:
    r"""
    Generate a list of random intervals.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: random_intervals(5) # random
        [(-oo, -4), [1, 2], (0, 1], (-oo, 1/2], (0, oo)]
    """
    return [
        random_interval(
            ring=ring,
            allow_infinity=allow_infinity,
            allow_open=allow_open,
            allow_empty=allow_empty,
        )
        for _ in range(length)
    ]


def intervals_from_sign_vector(sv) -> list[RealSet]:
    r"""
    Return intervals that are determined by a sign vector.

    EXAMPLES::

        sage: from sign_vectors import *
        sage: from vectors_in_intervals import *
        sage: intervals_from_sign_vector(sign_vector("++0-"))
        [(0, +oo), (0, +oo), {0}, (-oo, 0)]
    """
    lower_bounds = (
        0 if element > 0 else (-Infinity if element < 0 else 0) for element in sv
    )
    upper_bounds = (
        Infinity if element > 0 else (0 if element < 0 else 0) for element in sv
    )
    closed = [element == 0 for element in sv]
    return [
        interval_from_bounds(*bounds)
        for bounds in zip(lower_bounds, upper_bounds, closed, closed)
    ]

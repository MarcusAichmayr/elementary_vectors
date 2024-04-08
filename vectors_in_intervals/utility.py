"""Utility functions"""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.categories.sets_cat import EmptySetError
from sage.functions.other import floor, ceil
from sage.matrix.constructor import matrix
from sage.misc.prandom import randint
from sage.modules.free_module_element import vector
from sage.rings.continued_fraction import continued_fraction
from sage.rings.infinity import Infinity
from sage.rings.rational_field import QQ
from sage.sets.real_set import RealSet


def interval_from_bounds(
    lower_bound, upper_bound, lower_closed: bool = True, upper_closed: bool = True
) -> RealSet:
    r"""
    Construct an intervals.

    INPUT:

    - ``lower_bound`` -- a lower bound

    - ``upper_bound`` -- an upper bound

    - ``lower_closed`` -- a boolean (default: ``True``)

    - ``upper_closed`` -- a boolean (default: ``True``)

    OUTPUT:

    A ``RealSet`` object.

    - ``lower_bound`` and ``upper_bound`` are the left and right interval values.
      If ``lower_bound > upper_bound``, those elements will be exchanged.

    - ``lower_closed`` and ``upper_closed`` determine the intervals.

    - The left (or right) interval half of the interval is

        - closed if ``lower_closed`` (or ``upper_closed``) is ``True`` (default).

        - open if ``lower_closed`` (or ``upper_closed``) is ``False``.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import interval_from_bounds
        sage: interval_from_bounds(5, 6)
        [5, 6]
        sage: interval_from_bounds(6, 5, False, True)
        (5, 6]
        sage: interval_from_bounds(5, 5, False, True)
        {}
        sage: interval_from_bounds(-oo, 5)
        (-oo, 5]
        sage: interval_from_bounds(0, oo, False, False)
        (0, +oo)
    """
    if upper_bound < lower_bound:
        lower_bound, upper_bound = (upper_bound, lower_bound)
    if lower_bound == -Infinity:
        lower_closed = False
    if upper_bound == Infinity:
        upper_closed = False

    if lower_closed and upper_closed:
        interval = RealSet.closed(lower_bound, upper_bound)
    elif (not lower_closed) and (not upper_closed):
        interval = RealSet.open(lower_bound, upper_bound)
    elif lower_closed and (not upper_closed):
        interval = RealSet.closed_open(lower_bound, upper_bound)
    else:
        interval = RealSet.open_closed(lower_bound, upper_bound)

    return interval


def random_interval(
    ring=QQ,
    allow_infinity: bool = True,
    allow_open: bool = True,
    allow_empty: bool = False,
) -> RealSet:
    r"""Generate a random interval."""
    lower_bound = ring.random_element()
    upper_bound = ring.random_element()
    if allow_infinity:
        if randint(0, 5) == 0:
            lower_bound = -Infinity
        if randint(0, 5) == 0:
            upper_bound = Infinity
    if allow_open:
        lower_bound_closed = randint(0, 1) == 0
        upper_bound_closed = randint(0, 1) == 0
    else:
        lower_bound_closed = True
        upper_bound_closed = True
    interval = interval_from_bounds(
        lower_bound, upper_bound, lower_bound_closed, upper_bound_closed
    )
    if (not allow_empty) and interval.is_empty():
        return random_interval(ring, allow_infinity, allow_open, allow_empty)
    return interval


def simplest_element_in_interval(interval: RealSet):
    r"""
    Return the simplest rational element in an interval.

    INPUT:

    - ``interval`` -- an interval

    OUTPUT:
    If possible, an integer with smallest possible absolute value will be returned.
    Otherwise, a rational number with smallest possible denominator is constructed.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import simplest_element_in_interval
        sage: I = RealSet((1/2, oo))
        sage: simplest_element_in_interval(I)
        1
        sage: I = RealSet((-oo, 1/2))
        sage: simplest_element_in_interval(I)
        0
        sage: I = RealSet((-19, 0))
        sage: simplest_element_in_interval(I)
        -1
        sage: I = RealSet((0, 1))
        sage: simplest_element_in_interval(I)
        1/2
        sage: I = RealSet((-2/3, 0))
        sage: simplest_element_in_interval(I)
        -1/2
        sage: I = RealSet((4/3, 3/2))
        sage: simplest_element_in_interval(I)
        7/5
        sage: I = RealSet([0, 0])
        sage: simplest_element_in_interval(I)
        0
        sage: I = RealSet([5, 5])
        sage: simplest_element_in_interval(I)
        5
        sage: I = RealSet([sqrt(2), pi/2])
        sage: I
        [sqrt(2), 1/2*pi]
        sage: simplest_element_in_interval(I)
        3/2
        sage: I = RealSet([1/2, 1/2])
        sage: simplest_element_in_interval(I)
        1/2
    """
    if interval.is_empty():
        raise EmptySetError

    lower_bound = interval.inf()
    upper_bound = interval.sup()

    if (
        upper_bound - lower_bound > 1
        or floor(lower_bound) + 1 in interval
        or ceil(upper_bound) - 1 in interval
        or (
            lower_bound == upper_bound
            and (isinstance(lower_bound, int) or lower_bound.is_integer())
        )
    ):
        if 0 in interval:
            return 0
        if upper_bound == Infinity:
            if ceil(lower_bound) in interval:
                return ceil(lower_bound)
            return ceil(lower_bound) + 1
        if lower_bound == -Infinity:
            if floor(upper_bound) in interval:
                return floor(upper_bound)
            return floor(upper_bound) - 1
        if lower_bound == 0:
            return 1
        if upper_bound == 0:
            return -1
        if upper_bound < 0:
            return (
                floor(upper_bound)
                if floor(upper_bound) in interval
                else floor(upper_bound) - 1
            )
        # lower_bound > 0
        return (
            ceil(lower_bound)
            if ceil(lower_bound) in interval
            else ceil(lower_bound) + 1
        )
    return simplest_rational_in_interval(interval)


def simplest_rational_in_interval(interval: RealSet):
    r"""
    Find the rational with smallest denominator in a given interval.

    INPUT:

    - ``interval`` -- an interval that has no integer in it.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import simplest_rational_in_interval
        sage: I = RealSet((0, 1))
        sage: simplest_rational_in_interval(I)
        1/2
        sage: I = RealSet((1/3, 1))
        sage: simplest_rational_in_interval(I)
        1/2
        sage: I = RealSet((4/3, 2))
        sage: simplest_rational_in_interval(I)
        3/2
        sage: I = RealSet((1/3, 1/2))
        sage: simplest_rational_in_interval(I)
        2/5
        sage: I = RealSet((1/2, 241/287))
        sage: simplest_rational_in_interval(I)
        2/3
    """
    cfl = [
        floor(interval.inf()),
        2,
    ]  # continued fraction representation of inf(interval) + 1/2
    while True:
        value = continued_fraction(cfl).value()
        if value in interval:
            return value
        if value <= interval.inf():
            cfl = sb_child(cfl, left=False)
        else:
            cfl = sb_child(cfl, left=True)


def sb_child(cfl: list, left: bool):
    r"""
    Return a child of an element in the Stern-Brocot tree.

    INPUT:

    - ``cfl`` -- a list corresponding to a continued fraction

    - ``left`` -- a boolean

    OUTPUT:
    a list representing the continued fraction of a child of ``cfl``.
    If ``left`` is true, returns the left child of ``cfl``.
    Otherwise, the right child is returned.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import sb_child
        sage: sb_child([0, 2], True)
        [0, 3]
        sage: sb_child([0, 2], False)
        [0, 1, 2]
        sage: parent = continued_fraction_list(5/7)
        sage: parent
        [0, 1, 2, 2]
        sage: l_child = sb_child(parent, True)
        sage: l_child
        [0, 1, 2, 3]
        sage: continued_fraction(l_child).value()
        7/10
        sage: r_child = sb_child(parent, False)
        sage: r_child
        [0, 1, 2, 1, 2]
        sage: continued_fraction(r_child).value()
        8/11
    """
    if left == (len(cfl) % 2 == 0):
        return cfl[:-1] + [cfl[-1] + 1]
    return cfl[:-1] + [cfl[-1] - 1, 2]


def solve_left_for_roots(A, b):
    r"""
    Find a solution for ``x*A = b`` that works for matrices with roots.

    INPUT:

    - ``A`` -- a matrix

    - ``b`` -- a vector

    NOTE::

        The built in method ``solve_left`` for matrices does not always work.
    """
    M = matrix(list(A) + [-b]).T.right_kernel_matrix()
    x = matrix(M.column(-1)).solve_right(vector([1]))
    return (x * M)[:-1]

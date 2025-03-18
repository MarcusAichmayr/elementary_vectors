"""Interval classes."""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import annotations

from random import getrandbits

from sage.functions.other import floor, ceil
from sage.structure.sage_object import SageObject
from sage.rings.continued_fraction import continued_fraction
from sage.rings.infinity import minus_infinity
from sage.rings.infinity import Infinity
from sage.rings.rational_field import QQ

from sign_vectors import SignVector


class Interval(SageObject):
    r"""
    An interval.

    Also supports variables.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: Interval(0, 1)
        [0, 1]
        sage: Interval(0, 1, False, True)
        (0, 1]
        sage: Interval(0, 1, False, False)
        (0, 1)
        sage: Interval(0, 1, True, False)
        [0, 1)
        sage: Interval(0, 0)
        {0}
        sage: Interval(0, 0, False, False)
        {}
        sage: Interval(0, 0, True, True)
        {0}
        sage: Interval(0, 0, True, False)
        {}
        sage: Interval(0, 0, False, True)
        {}
        sage: Interval(-oo, 4)
        (-oo, 4]
        sage: Interval(-oo, 4, False, False)
        (-oo, 4)
        sage: I = Interval(-3, +oo, False, False)
        sage: I
        (-3, +oo)
        sage: 0 in I
        True
        sage: -3 in I
        False
        sage: Interval.random() # random
        (-1, 1)

    Variables are supported, too::

        sage: var('a')
        a
        sage: Interval(a, 1)
        [a, 1]
    """
    def __init__(self, lower, upper, lower_closed: bool = True, upper_closed: bool = True) -> None:
        if lower > upper:
            raise ValueError("The lower bound must be less than or equal to the upper bound.")

        if lower == upper and (not lower_closed or not upper_closed):
            lower = 0
            upper = 0
            lower_closed = False
            upper_closed = False

        self.lower = lower
        self.upper = upper
        if lower == minus_infinity:
            lower_closed = False
        if upper == Infinity:
            upper_closed = False
        self.lower_closed = lower_closed
        self.upper_closed = upper_closed

    def __contains__(self, x) -> bool:
        if self.lower_closed and x == self.lower:
            return True
        if self.upper_closed and x == self.upper:
            return True
        return self.lower < x < self.upper

    def _repr_(self) -> str:
        if self.is_empty():
            return "{}"
        if self.is_pointed():
            return f"{{{self.lower}}}"
        return f"{'[' if self.lower_closed else '('}{'-oo' if self.lower == minus_infinity else self.lower}, {'+oo' if self.upper == Infinity else self.upper}{']' if self.upper_closed else ')'}"

    def __eq__(self, other) -> bool:
        return self.lower == other.lower and self.upper == other.upper and self.lower_closed == other.lower_closed and self.upper_closed == other.upper_closed

    def is_empty(self) -> bool:
        r"""Return whether the interval is empty."""
        return self.lower == self.upper and (not self.lower_closed or not self.upper_closed)

    def is_open(self) -> bool:
        r"""Return whether the interval is open."""
        return self.lower_closed == self.upper_closed == False

    def is_closed(self) -> bool:
        r"""Return whether the interval is closed."""
        return self.lower_closed and self.upper_closed

    def is_unbounded(self) -> bool:
        r"""Return whether the interval is unbounded."""
        return self.lower == minus_infinity or self.upper == Infinity

    def is_bounded(self) -> bool:
        r"""Return whether the interval is bounded."""
        return not self.is_unbounded()

    def is_pointed(self) -> bool:
        r"""
        Return whether the interval is a point.
        
        EXAMPLES::
        
            sage: from vectors_in_intervals import *
            sage: Interval(0, 0).is_pointed()
            True
            sage: Interval(0, 1).is_pointed()
            False
        """
        return self.lower == self.upper and self.lower_closed and self.upper_closed

    def __bool__(self) -> bool:
        return not self.is_empty()

    def infimum(self):
        r"""
        Return the infimum of the interval.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval(0, 1).infimum()
            0
            sage: Interval(-oo, 0).infimum()
            -Infinity
            sage: Interval(0, 0, False, False).infimum()
            +Infinity
        """
        if self.is_empty():
            return Infinity
        return self.lower

    def supremum(self):
        r"""
        Return the supremum of the interval.
        
        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval(0, 1).supremum()
            1
            sage: Interval(0, +oo).supremum()
            +Infinity
            sage: Interval(0, 0, False, False).supremum()
            -Infinity
        """
        if self.is_empty():
            return minus_infinity
        return self.upper

    def an_element(self):
        r"""
        Return an element of the interval.

        EXAMPLES::
        
            sage: from vectors_in_intervals import *
            sage: Interval(0, 1).an_element()
            0
            sage: Interval(0, 1, False, True).an_element()
            1
            sage: Interval(0, 1, False, False).an_element()
            1/2
            sage: Interval(-oo, 0).an_element()
            0
            sage: Interval(-oo, +oo).an_element()
            0
            sage: Interval(5, +oo, False, False).an_element()
            6
            sage: Interval(0, 0, False, False).an_element()
            Traceback (most recent call last):
            ...
            ValueError: The interval is empty.
        """
        if self.is_empty():
            raise ValueError("The interval is empty.")
        if self.lower_closed:
            return self.lower
        if self.upper_closed:
            return self.upper
        if self.is_bounded():
            return (self.lower + self.upper) / 2
        if self.lower == minus_infinity:
            if self.upper == Infinity:
                return 0
            return self.upper - 1
        return self.lower + 1

    def simplest_element(self):
        r"""
        Return the simplest rational in this interval.

        OUTPUT:
        If possible, an integer with smallest possible absolute value will be returned.
        Otherwise, a rational number with smallest possible denominator is constructed.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval(1/2, +oo, False, False).simplest_element()
            1
            sage: Interval(-oo, 1/2, False, False).simplest_element()
            0
            sage: Interval(-19, 0, False, False).simplest_element()
            -1
            sage: Interval(0, 1, False, False).simplest_element()
            1/2
            sage: Interval(-2/3, 0, False, False).simplest_element()
            -1/2
            sage: Interval(4/3, 3/2, False, False).simplest_element()
            7/5
            sage: Interval(0, 0).simplest_element()
            0
            sage: Interval(5, 5).simplest_element()
            5
            sage: Interval(sqrt(2), pi/2).simplest_element()
            3/2
            sage: Interval(1/2, 1/2).simplest_element()
            1/2
        """
        if self.is_empty():
            raise ValueError("The interval is empty.")

        if self.is_pointed():
            return self.lower
        if 0 in self:
            return 0
        if (
            self.upper - self.lower > 1
            or floor(self.lower) + 1 in self
            or ceil(self.upper) - 1 in self
        ):
            if self.lower == minus_infinity:
                floor_upper = floor(self.upper)
                return floor_upper if floor_upper in self else floor_upper - 1
            if self.upper == Infinity:
                ceil_lower = ceil(self.lower)
                return ceil_lower if ceil_lower in self else ceil_lower + 1
            if self.lower == 0:
                return 1
            if self.upper == 0:
                return -1
            if self.upper < 0:
                floor_upper = floor(self.upper)
                return (
                    floor_upper if floor_upper in self else floor_upper - 1
                )
            # self.lower > 0
            ceil_lower = ceil(self.lower)
            return (
                ceil_lower if ceil_lower in self else ceil_lower + 1
            )
        return self._simplest_rational()

    def _simplest_rational(self):
        r"""
        Find the rational with smallest denominator in a given interval.

        INPUT:

        - ``interval`` -- an interval that has no integer in it.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval(0, 1, False, False)._simplest_rational()
            1/2
            sage: Interval(1/3, 1, False, False)._simplest_rational()
            1/2
            sage: Interval(4/3, 2, False, False)._simplest_rational()
            3/2
            sage: Interval(1/3, 1/2, False, False)._simplest_rational()
            2/5
            sage: Interval(1/2, 241/287, False, False)._simplest_rational()
            2/3
        """
        cfl = [floor(self.infimum()), 2]  # continued fraction representation of inf(interval) + 1/2

        def sb_child(cfl: list, left: bool):
            r"""
            Return a child of an element in the Stern-Brocot tree.

            INPUT:

            - ``cfl`` -- a list corresponding to a continued fraction

            - ``left`` -- a boolean
            """
            if left == (len(cfl) % 2 == 0):
                return cfl[:-1] + [cfl[-1] + 1]
            return cfl[:-1] + [cfl[-1] - 1, 2]

        while True:
            value = continued_fraction(cfl).value()
            if value in self:
                return value
            if value <= self.infimum():
                cfl = sb_child(cfl, left=False)
            else:
                cfl = sb_child(cfl, left=True)

    @staticmethod
    def random(ring=QQ, allow_infinity: bool = True, allow_empty: bool = False) -> Interval:
        r"""
        Generate a random interval.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval.random() # random
            (-1, 1)
        """
        if allow_infinity and getrandbits(3) == 0:
            lower = minus_infinity
        else:
            lower = ring.random_element()
        if allow_infinity and getrandbits(3) == 0:
            upper = Infinity
        else:
            upper = ring.random_element()
        if lower > upper:
            lower, upper = upper, lower
        lower_closed = bool(getrandbits(1))
        upper_closed = bool(getrandbits(1))

        interval = Interval(lower, upper, lower_closed, upper_closed)
        if not allow_empty and interval.is_empty():
            return Interval.random(ring, allow_infinity, allow_empty)
        return interval

    @staticmethod
    def open(lower, upper) -> Interval:
        r"""
        Return an open interval.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval.open(0, 1)
            (0, 1)
        """
        return Interval(lower, upper, False, False)

    @staticmethod
    def closed(lower, upper) -> Interval:
        r"""
        Return a closed interval.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval.closed(0, 1)
            [0, 1]
        """
        return Interval(lower, upper, True, True)

    @staticmethod
    def empty() -> Interval:
        r"""
        Return the empty interval.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Interval.empty()
            {}
        """
        return Interval(0, 0, False, False)


class Intervals(SageObject):
    r"""
    A Cartesian product of intervals.

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: Intervals([Interval(0, 1), Interval(-5, 2, False, False), Interval(0, 1)])
        [0, 1] x (-5, 2) x [0, 1]
        sage: vector([0, 1]) in Intervals([Interval(0, 1), Interval(-5, 2)])
        True
        sage: Intervals.random(3) # random
        [0, +oo) x (-5, 2) x (0, 1]
    """
    def __init__(self, intervals: list) -> None:
        self.intervals = intervals

    def __contains__(self, iterable) -> bool:
        return all(entry in interval for entry, interval in zip(iterable, self.intervals))

    def __len__(self) -> int:
        return len(self.intervals)

    def __getitem__(self, i) -> Interval:
        return self.intervals[i]

    def _repr_(self) -> str:
        if len(self) == 0:
            return "()"
        return " x ".join(str(interval) for interval in self.intervals)

    def __eq__(self, other) -> bool:
        return self.intervals == other.intervals

    def __hash__(self) -> int:
        return hash(tuple(self.intervals))

    def is_empty(self) -> bool:
        if len(self) == 0:
            return True
        return any(interval.is_empty() for interval in self.intervals)

    def is_open(self) -> bool:
        return all(interval.is_open() for interval in self.intervals)

    def is_closed(self) -> bool:
        return all(interval.is_closed() for interval in self.intervals)

    def is_half_open(self) -> bool:
        return any(interval.is_half_open() for interval in self.intervals)

    def is_unbounded(self) -> bool:
        return any(interval.is_unbounded() for interval in self.intervals)

    def is_bounded(self) -> bool:
        return not self.is_unbounded()

    def is_pointed(self) -> bool:
        return all(interval.is_pointed() for interval in self.intervals)

    def __bool__(self) -> bool:
        return not self.is_empty()

    def an_element(self):
        return [interval.an_element() for interval in self.intervals]

    def simplest_element(self):
        return [interval.simplest_element() for interval in self.intervals]

    def __iter__(self):
        return iter(self.intervals)

    @staticmethod
    def random(length: int, ring=QQ) -> Intervals:
        r"""
        Generate a random list of intervals.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Intervals.random(3) # random
            [0, +oo) x (-5, 2) x (0, 1]
        """
        return Intervals([Interval.random(ring) for _ in range(length)])

    @staticmethod
    def from_bounds(lower_bounds: list, upper_bounds: list, lower_bounds_closed: bool = True, upper_bounds_closed: bool = True) -> Intervals:
        r"""
        Return intervals that are determined by bounds.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: Intervals.from_bounds([0, -5, 0], [1, 2, +oo])
            [0, 1] x [-5, 2] x [0, +oo)
            sage: Intervals.from_bounds([0, -5, 0], [1, 2, +oo], False, False)
            (0, 1) x (-5, 2) x (0, +oo)
            sage: Intervals.from_bounds([0, -5, 0], [1, 2, +oo], True, False)
            [0, 1) x [-5, 2) x [0, +oo)
            sage: Intervals.from_bounds([0, -5, 0], [1, 2, +oo], [True, False, True], [True, True, False])
            [0, 1] x (-5, 2] x [0, +oo)
        """
        length = len(lower_bounds)
        if lower_bounds_closed is True:
            lower_bounds_closed = [True] * length
        elif lower_bounds_closed is False:
            lower_bounds_closed = [False] * length

        if upper_bounds_closed is True:
            upper_bounds_closed = [True] * length
        elif upper_bounds_closed is False:
            upper_bounds_closed = [False] * length

        return Intervals([
            Interval(*bounds)
            for bounds in zip(
                lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed
            )
        ])

    @staticmethod
    def from_sign_vector(sv: SignVector) -> Intervals:
        r"""
        Return intervals that are determined by a sign vector.

        EXAMPLES::

            sage: from vectors_in_intervals import *
            sage: from sign_vectors import *
            sage: sv = sign_vector("+0-")
            sage: Intervals.from_sign_vector(sv)
            (0, +oo) x {0} x (-oo, 0)
        """
        return Intervals([
            Interval(
                0 if element > 0 else (minus_infinity if element < 0 else 0),
                Infinity if element > 0 else (0 if element < 0 else 0),
                element == 0,
                element == 0
            )
            for element in sv
        ])

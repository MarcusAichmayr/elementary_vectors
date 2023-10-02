r"""
Finding vectors with components in intervals.

With this module, we can check whether there is a vector in a subspace such that
the components lie in given intervals.

There is also an algorithmic approach to construct such a vector.

EXAMPLES:

First, we load the functions from the package::

    sage: from vectors_in_intervals import *

The package offers the function :func:`vectors_in_intervals.setup_intervals.intervals_from_bounds`
that helps us creating lists of intervals.
This function takes two lists as input,
the corresponding elements in those lists determine the intervals::

    sage: I = intervals_from_bounds([0, 1, -1], [1, 2, -1])
    sage: I
    [[0, 1], [1, 2], {-1}]

Next, we define a vector::

    sage: v = vector([1,1,1])

Is there a normal vector, such that the components lie in the intervals defined above?
We call :func:`vectors_in_intervals.existence.exists_orthogonal_vector` to answer this question::

    sage: exists_orthogonal_vector(v, I)
    True

This vector can be constructed by
:func:`vectors_in_intervals.construction.construct_orthogonal_vector`::

    sage: construct_orthogonal_vector(v, I)
    (0, 1, -1)

We define another vector. This time, there is no solution::

    sage: v = vector([1,1,-1])
    sage: exists_orthogonal_vector(v, I)
    False
    sage: construct_orthogonal_vector(v, I)
    Traceback (most recent call last):
    ...
    ValueError: There is no solution.

Next, we consider open intervals::

    sage: I = intervals_from_bounds([0, 1, -1], [1, 2, 1], False, [True, True, False])
    sage: I
    [(0, 1], (1, 2], (-1, 1)]
    sage: v = vector([1,0,1])
    sage: exists_orthogonal_vector(v, I)
    True
    sage: construct_orthogonal_vector(v, I)
    (1/2, 2, -1/2)

We can even consider unbounded intervals::

    sage: I = intervals_from_bounds([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
    sage: I
    [(0, +oo), (1, 2], (-oo, -2)]
    sage: v = vector([-1,1,-1])
    sage: exists_orthogonal_vector(v, I)
    True
    sage: construct_orthogonal_vector(v, I)
    (5, 2, -3)

The most important functions of this module are :func:`vectors_in_intervals.existence.exists_vector`
and :func:`vectors_in_intervals.construction.construct_vector`.
Given a matrix ``M`` and a list of intervals,
we want to examine whether there exists a vector in the row space of ``M``,
such that the components lie in the given intervals::

    sage: M = matrix([1, 1, 0])
    sage: lower_bounds = [2, 5, -1]
    sage: upper_bounds = [5, 6, 1]

First, we consider closed intervals::

    sage: I = intervals_from_bounds(lower_bounds, upper_bounds)
    sage: I
    [[2, 5], [5, 6], [-1, 1]]
    sage: exists_vector(M, I)
    True
    sage: construct_vector(M, I)
    (5, 5, 0)

Next, we take open intervals. This time, there is no solution::

    sage: I = intervals_from_bounds(lower_bounds, upper_bounds, False, False)
    sage: I
    [(2, 5), (5, 6), (-1, 1)]
    sage: exists_vector(M, I)
    False
    sage: construct_vector(M, I)
    Traceback (most recent call last):
    ...
    ValueError: There is no solution.

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
    sage: construct_vector(M, I)
    (2, 5, 7, 5)
"""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import absolute_import

from .setup_intervals import is_vector_in_intervals, intervals_from_bounds, intervals_from_sign_vector
from .existence import exists_vector, exists_orthogonal_vector
from .construction import construct_vector, construct_orthogonal_vector, vector_from_sign_vector, sign_vectors_in_intervals

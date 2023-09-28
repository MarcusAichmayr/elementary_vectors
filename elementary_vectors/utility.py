"""Utility functions."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.categories.sets_cat import EmptySetError
from sage.functions.other import floor, ceil
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.continued_fraction import continued_fraction
from sage.rings.infinity import Infinity
from sage.sets.real_set import RealSet
from sage.symbolic.ring import SR

from sign_vectors import sign_vector


def sign_determined(expression):
    r"""
    Check whether the sign of a number or symbolic expression ``expression`` is uniquely determined.

    EXAMPLES::

        sage: from elementary_vectors.utility import sign_determined

    Integers have always a unique sign::

        sage: sign_determined(2)
        True
        sage: sign_determined(-5)
        True

    Now, we consider a variable::

        sage: var('a')
        a
        sage: sign_determined(a)
        False
        sage: assume(a >= 0)
        sage: sign_determined(a)
        False
        sage: assume(a != 0)
        sage: sign_determined(a)
        True
        sage: sign_determined(a - 1)
        False
    """
    return bool(SR(expression) > 0 or SR(expression) < 0 or SR(expression) == 0)


def kernel_vector_support_given(M, indices):
    r"""
    Return a vector in the right kernel of ``M`` such that the support is a subset of ``indices``.

    INPUT:

    - ``M`` -- a matrix

    - ``indices`` -- a list of indices

    OUTPUT:
    a vector ``v`` in the right kernel of ``M`` such that
    the support is a subset of ``indices``.

    EXAMPLES::

        sage: from elementary_vectors.utility import kernel_vector_support_given
        sage: M = matrix([[1, 2, 0, 0], [0, 1, -1, 0]])
        sage: kernel_vector_support_given(M, [0, 1, 2])
        (2, -1, -1, 0)
        sage: kernel_vector_support_given(M, [3])
        (0, 0, 0, 1)
        sage: kernel_vector_support_given(M, [0, 3])
        (0, 0, 0, 1)
    """
    M_I = M.matrix_from_columns(indices)
    try:
        subvector_list = list(M_I.right_kernel_matrix()[0])
    except IndexError:
        raise ValueError("Right kernel of ``M`` restricted to the columns ``indices`` is empty.")
    for k in range(M.ncols()):
        if not k in indices:
            subvector_list.insert(k, 0)
    return vector(subvector_list)


def interval_from_bounds(lower_bound, upper_bound, lower_closed=True, upper_closed=True):
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

        sage: from elementary_vectors.vectors_in_intervals import interval_from_bounds
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


def simplest_element_in_interval(interval):
    r"""
    Return the simplest rational element in an interval.

    INPUT:

    - ``interval`` -- an interval (``RealSet``)

    OUTPUT:
    If possible, an integer with smallest possible absolute value will be returned.
    Otherwise, a rational number with smallest possible denominator is constructed.

    EXAMPLES::

        sage: from elementary_vectors.utility import simplest_element_in_interval
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
        or (lower_bound.is_integer() and lower_bound == upper_bound)
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
            return floor(upper_bound) if floor(upper_bound) in interval else floor(upper_bound) - 1
        # lower_bound > 0
        return ceil(lower_bound) if ceil(lower_bound) in interval else ceil(lower_bound) + 1
    return simplest_rational_in_interval(interval)


def simplest_rational_in_interval(interval):
    r"""
    Find the rational with smallest denominator in a given interval.

    INPUT:

    - ``interval`` -- an interval (``RealSet``) that has no integer in it.

    EXAMPLES::

        sage: from elementary_vectors.utility import simplest_rational_in_interval
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
    cfl = [floor(interval.inf()), 2] # continued fraction representation of inf(interval) + 1/2
    while True:
        value = continued_fraction(cfl).value()
        if value in interval:
            return value
        if value <= interval.inf():
            cfl = sb_child(cfl, left=False)
        else:
            cfl = sb_child(cfl, left=True)


def sb_child(cfl, left):
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

        sage: from elementary_vectors.utility import sb_child
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


def solve_left(A, b):
    r"""
    Find a solution for ``x*A = b`` that works for matrizes with roots.
    
    INPUT:
    
    - ``A`` -- a matrix
    
    - ``b`` -- a vector
    
    NOTE::
    
        The built in method ``solve_left`` for matrices does not always work.
    """
    M = matrix(list(A) + [-b]).T.right_kernel_matrix()
    x = matrix(M.column(-1)).solve_right(vector([1]))
    return (x * M)[:-1]


def conformal_elimination(vector1, vector2, indices=None):
    r"""
    Apply conformal elimination to two real vectors to find a new vector.

    INPUT:

    - ``vector1`` -- a real vector

    - ``vector2`` -- a real vector

    - ``indices`` -- a list of indices (default: ``[]``)

    OUTPUT:

    Returns a new vector ``z = x + a y`` where ``a > 0``, such that ``z[e] == 0``
    for some ``e`` in ``indices`` and ``Z_k <= X_k`` for ``k`` in ``indices``
    and ``Z_f = (X o Y)_f`` for ``f`` not in ``D(X, Y)``.
    Here, ``X``, ``Y`` and ``Z`` are the sign vectors
    corresponding to ``x``, ``y`` and ``z``.

    .. NOTE::

        If ``indices`` is not given, the whole list of separating elements
        will be considered instead. (default)
    
    EXAMPLES::
    
        sage: from elementary_vectors.utility import conformal_elimination
        sage: x = vector([1, 0, 2])
        sage: y = vector([-1, 1, 1])
        sage: conformal_elimination(x, y)
        (0, 1, 3)
    """
    if indices is None:
        indices = []
    if vector1.length() != vector2.length():
        raise ValueError('Vectors have different length.')
    separating_elements = sign_vector(vector1).separating_elements(sign_vector(vector2))
    if separating_elements == []:
        raise ValueError('List of separating elements is empty.')
    if indices == []:
        indices = separating_elements
    elif not all(s in separating_elements for s in indices):
        raise ValueError('indices is not a subset of separating_elements.')
    lam = max(vector1[e] / vector2[e] for e in indices)  # lam < 0 since e in separating_elements
    return vector1 - lam * vector2

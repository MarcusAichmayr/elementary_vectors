"""Utility functions."""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.symbolic.ring import SR
from sign_vectors import sign_vector
from sage.modules.free_module_element import vector
from sage.sets.real_set import RealSet
from sage.rings.infinity import Infinity
from sage.functions.other import floor, ceil
from sage.rings.continued_fraction import continued_fraction


def sign_determined(a):
    r"""
    Check whether the sign of a number or symbolic expression ``a`` is uniquely determined.

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
    return bool(SR(a) > 0 or SR(a) < 0 or SR(a) == 0)


def vector_from_matrix(M, I):
    r"""
    Return a vector in the right kernel of ``M`` such that
    the support is a subset of ``I``.

    INPUT:

    - ``M`` -- a matrix

    - ``I`` -- a list of indices

    OUTPUT:
    a vector ``v`` in the right kernel of ``M`` such that
    the support is a subset of ``I``.

    EXAMPLES::

        sage: from elementary_vectors.utility import vector_from_matrix
        sage: M = matrix([[1, 2, 0, 0], [0, 1, -1, 0]])
        sage: vector_from_matrix(M, [0, 1, 2])
        (2, -1, -1, 0)
        sage: vector_from_matrix(M, [3])
        (0, 0, 0, 1)
        sage: vector_from_matrix(M, [0, 3])
        (0, 0, 0, 1)
    """
    n = M.ncols()
    M_I = M.matrix_from_columns(I)
    try:
        l = list(M_I.right_kernel_matrix()[0])
    except IndexError:
        raise(ValueError("Right kernel of ``M`` restricted to the columns ``I`` is empty."))
    for k in range(n):
        if not k in I:
            l.insert(k, 0)
    return vector(l)


def setup_interval(L, R, l=True, r=True):
    r"""
    Construct an intervals.

    INPUT:

    - ``L`` -- a lower bound

    - ``R`` -- an upper bound

    - ``l`` -- a boolean (default: ``True``)

    - ``r`` -- a boolean (default: ``True``)

    OUTPUT:

    A ``RealSet`` objects.

    - ``L`` and ``R`` are the left and right interval values.
      If ``L > R``, those elements will be exchanged.

    - ``l`` and ``r`` determine the intervals.

    - The left (or right) interval half of the interval is

        - closed if ``l`` (or ``r``) is ``True`` (default).

        - open if ``l`` (or ``r``) is ``False``.

    EXAMPLES::

        sage: from elementary_vectors.vectors_in_intervals import setup_interval
        sage: setup_interval(5, 6)
        [5, 6]
        sage: setup_interval(6, 5, False, True)
        (5, 6]
        sage: setup_interval(5, 5, False, True)
        {}
        sage: setup_interval(-oo, 5)
        (-oo, 5]
        sage: setup_interval(0, oo, False, False)
        (0, +oo)
    """
    if R < L:
        L, R = (R, L)
    if L == -Infinity:
        l = False
    if R == Infinity:
        r = False

    if l and r:
        interval = RealSet.closed(L, R)
    elif (not l) and (not r):
        interval = RealSet.open(L, R)
    elif l and (not r):
        interval = RealSet.closed_open(L, R)
    else:
        interval = RealSet.open_closed(L, R)

    return interval


def simplest_element_in_interval(I):
    r"""
    Return the simplest rational element in an interval.

    INPUT:

    - ``I`` -- an interval (``RealSet``)

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
    if I.is_empty():
        raise EmptySetError

    a = I.inf()
    b = I.sup()
    if (
        b - a > 1
        or floor(a) + 1 in I
        or ceil(b) - 1 in I
        or (a.is_integer() and a == b)
    ):
        if 0 in I:
            return 0
        elif b == Infinity:
            if ceil(a) in I:
                return ceil(a)
            else:
                return ceil(a) + 1
        elif a == -Infinity:
            if floor(b) in I:
                return floor(b)
            else:
                return floor(b) - 1
        else:
            if a == 0:
                return 1
            elif b == 0:
                return -1
            elif b < 0:
                return floor(b) if floor(b) in I else floor(b) - 1
            else: # a > 0
                return ceil(a) if ceil(a) in I else ceil(a) + 1
    else:
        return simplest_rational_in_interval(I)


def simplest_rational_in_interval(I):
    r"""
    Find the rational with smallest denominator in a given interval.

    INPUT:

    - ``I`` -- an interval (``RealSet``) that has no integer in it.

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
    cfl = [floor(I.inf()), 2] # continued fraction representation of inf(I) + 1/2
    while True:
        val = continued_fraction(cfl).value()
        if val in I:
            return val
        else:
            if val <= I.inf():
                cfl = sb_child(cfl, left=False)
            else: # >= I.sup()
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
    else:
        return cfl[:-1] + [cfl[-1] - 1, 2]


def conformal_elimination(x, y, S=None):
    r"""
    Apply conformal elimination to two real vectors to find a new vector.

    INPUT:

    - ``x`` -- a real vector

    - ``y`` -- a real vector

    - ``S`` -- a list of indices (default: ``[]``)

    OUTPUT:

    Returns a new vector ``z = x + a y`` where ``a > 0``, such that ``z[e] == 0``
    for some ``e`` in ``S`` and ``Z_S <= X_S`` and ``Z_f = (X o Y)_f`` for ``f``
    not in ``D(X, Y)``. Here, ``X``, ``Y`` and ``Z`` are the sign vectors
    corresponding to ``x``, ``y`` and ``z``.

    .. NOTE::

        If ``S`` is the empty list ``[]``, the whole list of separating elements
        will be considered instead. (default)
    """
    if S is None:
        S = []
    if x.length() != y.length():
        raise ValueError('Vectors have different length.')
    X = sign_vector(x)
    Y = sign_vector(y)
    D = X.separating_elements(Y)
    if D == []:
        raise ValueError('List of separating elements is empty.')
    if S == []:
        S = D
    elif not all(s in D for s in S):
        raise ValueError('S is not a subset of D.')
    lam = max([x[e]/y[e] for e in S])  # x[e]/y[e] < 0 since e in D
    return x - lam*y

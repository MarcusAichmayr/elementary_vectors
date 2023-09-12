r"""
Finding vectors in intervals.

With this module, we can check whether there is a vector in a subspace such that
the components lie in given intervals.

There is also an algorithmic approach to construct such a vector.

EXAMPLES:

First, we load the functions from the package::

    sage: from elementary_vectors import *

The package offers the function :func:`~setup_intervals` that helps us creating lists of intervals.
This function takes two lists as input,
the corresponding elements in those lists determine the intervals::

    sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
    sage: I
    [[0, 1], [1, 2], {-1}]

Next, we define a vector::

    sage: v = vector([1,1,1])

Is there a normal vector, such that the components lie in the intervals defined above?
We call :func:`~exists_orthogonal_vector` to answer this question::

    sage: exists_orthogonal_vector(v, I)
    True

This vector can be constructed by :func:`~construct_normal_vector`::

    sage: construct_normal_vector(v, I)
    (0, 1, -1)

We define another vector. This time, there is no solution::

    sage: v = vector([1,1,-1])
    sage: exists_orthogonal_vector(v, I)
    False
    sage: construct_normal_vector(v, I)
    Traceback (most recent call last):
    ...
    ValueError: There is no solution.

Next, we consider open intervals::

    sage: I = setup_intervals([0, 1, -1], [1, 2, 1], False, [True, True, False])
    sage: I
    [(0, 1], (1, 2], (-1, 1)]
    sage: v = vector([1,0,1])
    sage: exists_orthogonal_vector(v, I)
    True
    sage: construct_normal_vector(v, I)
    (1/2, 2, -1/2)

We can even consider unbounded intervals::

    sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
    sage: I
    [(0, +oo), (1, 2], (-oo, -2)]
    sage: v = vector([-1,1,-1])
    sage: exists_orthogonal_vector(v, I)
    True
    sage: construct_normal_vector(v, I)
    (5, 2, -3)

The most important functions of this module are
:func:`~exists_vector` and :func:`~construct_vector`.
Given a matrix ``M`` and a list of intervals,
we want to examine whether there exists a vector in the rowspace of ``M``,
such that the components lie in the given intervals::

    sage: M = matrix([1, 1, 0])
    sage: lower_bounds = [2, 5, -1]
    sage: upper_bounds = [5, 6, 1]

First, we consider closed intervals::

    sage: I = setup_intervals(lower_bounds, upper_bounds)
    sage: I
    [[2, 5], [5, 6], [-1, 1]]
    sage: exists_vector(M, I)
    True
    sage: construct_vector(M, I)
    (5, 5, 0)

Next, we take open intervals. This time, there is no solution::

    sage: I = setup_intervals(lower_bounds, upper_bounds, False, False)
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
    sage: I = setup_intervals(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
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

from .functions import elementary_vectors
from .utility import simplest_element_in_interval, setup_interval, solve_left
from sage.sets.real_set import RealSet
from sage.rings.infinity import Infinity
from sage.calculus.var import var
from sage.symbolic.relation import solve
from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.constructor import matrix


def setup_intervals(lower_bounds, upper_bounds, lower_bounds_closed=True, upper_bounds_closed=True):
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

        sage: from elementary_vectors.vectors_in_intervals import setup_intervals
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    By default, the intervals are closed::

        sage: setup_intervals(lower_bounds, upper_bounds)
        [[2, 5], [5, 6], [-1, 1]]

    We obtain open intervals if both ``lower_bounds_closed`` and ``upper_bounds_closed`` are false::

        sage: setup_intervals(lower_bounds, upper_bounds, False, False)
        [(2, 5), (5, 6), (-1, 1)]

    Mixed intervals are also possible::

        sage: setup_intervals(lower_bounds, upper_bounds, False)
        [(2, 5], (5, 6], (-1, 1]]
        sage: setup_intervals(lower_bounds, upper_bounds, [False, True, True], [True, True, False])
        [(2, 5], [5, 6], [-1, 1)]

    We can also specify unbounded intervals.
    Note that bounds at infinity are always open::

        sage: lower_bounds = [-oo, 2, -oo]
        sage: upper_bounds = [-5, oo, oo]
        sage: setup_intervals(lower_bounds, upper_bounds)
        [(-oo, -5], [2, +oo), (-oo, +oo)]

    Finite and empty intervals are represented as usual::

        sage: lower_bounds = [0, 2]
        sage: upper_bounds = [0, 2]
        sage: setup_intervals(lower_bounds, upper_bounds, [True, False], True)
        [{0}, {}]

    If the lower bound is greater than the upper bound, those bounds will be exchanged::

        sage: lower_bounds = [1, 4]
        sage: upper_bounds = [2, 3]
        sage: setup_intervals(lower_bounds, upper_bounds, False)
        [(1, 2], (3, 4]]
    """
    n = len(lower_bounds)
    if len(upper_bounds) != n:
        raise ValueError('``upper_bounds`` should be a list of length ' + str(n) + '.')

    if lower_bounds_closed is True:
        lower_bounds_closed = [True] * n
    elif lower_bounds_closed is False:
        lower_bounds_closed = [False] * n
    elif len(lower_bounds_closed) != n:
        raise ValueError('``lower_bounds_closed`` should be a list of length ' + str(n) + '.')

    if upper_bounds_closed is True:
        upper_bounds_closed = [True] * n
    elif upper_bounds_closed is False:
        upper_bounds_closed = [False] * n
    elif len(upper_bounds_closed) != n:
        raise ValueError('``upper_bounds_closed`` should be a list of length ' + str(n) + '.')

    return [setup_interval(*bounds) for bounds in zip(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)]


def simplest_vector_in_intervals(intervals):
    r"""
    Return the simplest vector such that each component is in a given interval.

    INPUT:

    - ``intervals`` -- a list of intervals (``RealSet``)

    OUTPUT:
    A vector with components in the intervals.
    If possible, each component is the integer
    with smallest possible absolute value in the corresponding interval.
    Otherwise, components are rational with smallest possible denominator.

    EXAMPLES::

        sage: from elementary_vectors.vectors_in_intervals import *
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]
        sage: I = setup_intervals(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: simplest_vector_in_intervals(I)
        (2, 5, 0)
        sage: I = setup_intervals(lower_bounds, upper_bounds, False, False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: simplest_vector_in_intervals(I)
        (3, 11/2, 0)
    """
    return vector(simplest_element_in_interval(interval) for interval in intervals)


def multiple_in_intervals_candidates(v, intervals):
    r"""
    Return the biggest interval ``J`` such that ``a v`` lies in given intervals for all ``a`` in ``J``.

    INPUT:

    - ``v`` -- a vector

    - ``intervals`` -- a list of intervals (``RealSet``)

    .. SEEALSO::

        :func:`~multiple_in_intervals`

    EXAMPLES::

        sage: from elementary_vectors.vectors_in_intervals import *
        sage: I = setup_intervals([0, 0, -1/4, -oo], [2, oo, 1/5, 1/9], False, [True, False, False, True])
        sage: v = vector([1, 5, -2, 0])
        sage: multiple_in_intervals_candidates(v, I)
        (0, 1/8)
        sage: v = vector([0, 5, -2, 0])
        sage: multiple_in_intervals_candidates(v, I)
        {}
        sage: I = setup_intervals([-2, -1/4, -7], [0, 1, 0], [True, False, False], [False, False, True])
        sage: v = vector([-2, -2, -1])
        sage: multiple_in_intervals_candidates(v, I)
        (0, 1/8)
        sage: v = vector([0, -2, 1])
        sage: I = setup_intervals([-2, 0, -1], [1, 1/3, 0], True, False)
        sage: multiple_in_intervals_candidates(v, I)
        (-1/6, 0)
    """
    a_inf = -Infinity
    a_sup = Infinity
    lower_bound_closed = False
    upper_bound_closed = False

    for entry, interval in zip(v, intervals):
        if entry == 0:
            if not 0 in interval:
                return RealSet()
        else:
            val_inf = interval.inf() / entry
            val_sup = interval.sup() / entry
            if entry > 0:
                if val_inf > a_inf:
                    a_inf = val_inf
                    lower_bound_closed = interval.inf() in interval
                elif val_inf == a_inf:
                    lower_bound_closed = interval.inf() in interval and lower_bound_closed
                if val_sup < a_sup:
                    a_sup = val_sup
                    upper_bound_closed = interval.sup() in interval
                elif val_sup == a_sup:
                    upper_bound_closed = interval.sup() in interval and upper_bound_closed
            else: # entry < 0
                if val_sup > a_inf:
                    a_inf = val_sup
                    lower_bound_closed = interval.sup() in interval
                elif val_sup == a_inf:
                    lower_bound_closed = interval.sup() in interval and lower_bound_closed
                if val_inf < a_sup:
                    a_sup = val_inf
                    upper_bound_closed = interval.inf() in interval
                elif val_inf == a_sup:
                    upper_bound_closed = interval.inf() in interval and upper_bound_closed
            if a_inf > a_sup:
                return RealSet()

    return setup_interval(a_inf, a_sup, lower_bound_closed, upper_bound_closed)


def multiple_in_intervals(v, intervals):
    r"""
    Return a multiple of a vector that lies in given intervals if possible.

    INPUT:

    - ``v`` -- a vector

    - ``intervals`` -- a list of intervals (``RealSet``)

    OUTPUT:
    Computes a multiple ``a v`` that lies in the intervals.
    The number ``a`` is chosen as simple as possible.
    However, there might exist a simpler multiple lying in the intervals.

    .. SEEALSO::

        :func:`.utility.simplest_element_in_interval`

    EXAMPLES::

        sage: from elementary_vectors.vectors_in_intervals import *
        sage: I = setup_intervals([0, 0, -1/4, -oo], [2, oo, 1/5, 1/9], False, [True, False, False, True])
        sage: v = vector([1, 5, -2, 0])
        sage: multiple_in_intervals(v, I)
        (1/9, 5/9, -2/9, 0)
    """
    return simplest_element_in_interval(multiple_in_intervals_candidates(v, intervals)) * v


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

    .. SEEALSO::

        :func:`~setup_intervals`

    EXAMPLES:

    We define several lists of intervals and vectors
    and apply the function::

        sage: from elementary_vectors import exists_orthogonal_vector, setup_intervals
        sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1, 1, 1])
        sage: exists_orthogonal_vector(v, I)
        True
        sage: v = vector([1, 1, -1])
        sage: exists_orthogonal_vector(v, I)
        False

    Next, we consider open intervals::

        sage: I = setup_intervals([0, 1, -1], [1, 2, 1], False, [True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1, 0, 1])
        sage: exists_orthogonal_vector(v, I)
        True

    We can even consider unbounded intervals::

        sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1, 1, -1])
        sage: exists_orthogonal_vector(v, I)
        True

    TESTS::

        sage: I = setup_intervals([-oo, 0], [oo, 0], False, False)
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


def exists_vector(data, intervals, certificate=False):
    r"""
    Return whether a vector exists in a given vector space such that the components lie in the specified intervals.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                  elementary vectors of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)
    
    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:

    Return whether there exists a vector in the vector space determined by data
    such that the components lie in the specified intervals.

    - If ``data`` is a matrix, then the elementary vectors in the kernel of ``M`` are computed.

    - If ``data`` is a list of elementary vectors, then those will be used.

    - If ``certificate`` is true and no vector exists, an elementary vector certifying non-existence is returned.

    ALGORITHM:

    The underlying algorithm is based on Minty's Lemma. (see [Min74]_)

    .. [Min74] Minty, G. J.:
       „A `from scratch` proof of a theorem of Rockafellar and Fulkerson“.
       In: Mathematical Programming 7 (1974), pp. 368-375.

    .. SEEALSO::

        :func:`~setup_intervals`
        :func:`~exists_orthogonal_vector`

    EXAMPLES::

        sage: from elementary_vectors import exists_vector, setup_intervals
        sage: M = matrix([1, 1, 0])
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    First, we consider closed intervals::

        sage: I = setup_intervals(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: exists_vector(M, I)
        True

    Next, we take open intervals::

        sage: I = setup_intervals(lower_bounds, upper_bounds, False, False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: exists_vector(M, I)
        False

    Since no vector exists, there is an elementary vector certifying this.
    To find one, we pass ``certificate=True``::
    
        sage: exists_vector(M, I, certificate=True)
        (1, -1, 0)

    Mixed intervals are also possible::

        sage: lower_bounds_closed = [True, True, False]
        sage: upper_bounds_closed = [False, True, True]
        sage: I = setup_intervals(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
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
        sage: I = setup_intervals(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: exists_vector(M, I)
        True

    TESTS::

        sage: I = setup_intervals([0, 0], [1, 0], False)
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
            if certificate:
                return v
            else:
                return False
    return True


def construct_normal_vector(v, intervals):
    r"""
    Construct a vector that is normal on a given vector and lies in the specified intervals.

    INPUT:

    - ``v`` -- a vector of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    OUTPUT:

    Returns a (rational) vector ``z`` such that the scalar product of ``z`` and ``v`` is zero
    and each component of ``z`` lies in the respective interval of the list ``intervals``.
    If no such vector exists, raises a ``ValueError`` instead.

    .. SEEALSO::

        :func:`~setup_intervals`
        :func:`~exists_orthogonal_vector`

    EXAMPLES::

        sage: from elementary_vectors import construct_normal_vector, setup_intervals
        sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1, 1, 1])
        sage: construct_normal_vector(v, I)
        (0, 1, -1)

    We define another vector. This time, there is no solution::

        sage: v = vector([1, 1, -1])
        sage: construct_normal_vector(v, I)
        Traceback (most recent call last):
        ...
        ValueError: There is no solution.

    Next, we consider open intervals::

        sage: I = setup_intervals([0, 1, -1], [1, 2, 1], False, [True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1, 0, 1])
        sage: construct_normal_vector(v, I)
        (1/2, 2, -1/2)

    We can even consider unbounded intervals::

        sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1, 1, -1])
        sage: construct_normal_vector(v, I)
        (5, 2, -3)
    """
    if not exists_orthogonal_vector(v, intervals):
        raise ValueError("There is no solution.")

    # construct z_min and z_max
    z_min = []
    z_max = []
    var('eps')
    var('lam')

    for entry, interval in zip(v, intervals):
        if entry == 0:
            z_min.append(simplest_element_in_interval(interval))
            z_max.append(simplest_element_in_interval(interval))
        else:
            lower_bound_closed = (interval.inf() + (0 if interval.inf() in interval else eps)) if interval.inf() != -Infinity else -lam
            upper_bound_closed = (interval.sup() + (0 if interval.sup() in interval else -eps)) if interval.sup() != Infinity else lam

            if entry > 0:
                z_min.append(lower_bound_closed)
                z_max.append(upper_bound_closed)
            else:
                z_min.append(upper_bound_closed)
                z_max.append(lower_bound_closed)

    z_min = vector(z_min)
    z_max = vector(z_max)

    # find candidates for eps and lam
    eps_candidates = []
    lam_candidates = []
    for entry, interval in zip(v, intervals):
        if entry != 0:
            if interval.inf() == -Infinity and interval.sup() == Infinity:  # (-oo, oo)
                lam_candidates.append(0)
            elif interval.inf() == -Infinity:
                if interval.sup() in interval:  # (-oo, b]
                    lam_candidates.append(-interval.sup())
                else:  # (-oo, b)
                    lam_candidates.append(1 - interval.sup())
                    eps_candidates.append(1)
            elif interval.sup() == Infinity:
                if interval.inf() in interval:  # [a, oo)
                    lam_candidates.append(interval.inf())
                else:  # (a, oo)
                    lam_candidates.append(1 + interval.inf())
                    eps_candidates.append(1)
            elif not interval.is_closed():
                if interval.sup() in interval or interval.inf() in interval:  # [a, b) or (a, b]
                    eps_candidates.append(interval.sup() - interval.inf())
                else:  # (a, b)
                    eps_candidates.append((interval.sup() - interval.inf())/2)

    if eps_candidates:
        for product in [z_min * v, z_max * v]:
            try:
                if eps in product.variables() and lam not in product.variables():
                    eps_candidates.append(solve(product, eps, solution_dict=True)[0][eps])
            except AttributeError:
                pass
        eps_min = min(eps_candidates)
        try:
            z_min = z_min(eps=eps_min)
        except TypeError:
            pass
        try:
            z_max = z_max(eps=eps_min)
        except TypeError:
            pass

    if lam_candidates:
        for product in [z_min * v, z_max * v]:
            if product != 0:
                try:
                    lam_candidates.append(solve(product, lam, solution_dict=True)[0][lam])
                except (IndexError, TypeError):
                    pass
        lam_max = max(lam_candidates)
        try:
            z_min = z_min(lam=lam_max)
        except TypeError:
            pass
        try:
            z_max = z_max(lam=lam_max)
        except TypeError:
            pass

    product_min = v * z_min
    product_max = v * z_max
    if product_min == 0:
        return z_min
    if product_max == 0:
        return z_max

    return (product_max * z_min - product_min * z_max) / (product_max - product_min)


def construct_vector(M, intervals, evs=None):
    r"""
    Return a vector of a given vectorspace such that the components lie in given intervals.

    INPUT:

    - ``M`` -- a matrix with ``n`` columns

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    - ``evs`` -- an optional iterable of elementary vectors

    OUTPUT:

    Returns a vector in the rowspace of ``M`` such that each component lies
    in the respective interval of the list ``intervals``.
    If no such vector exists, raises a ``ValueError`` instead.

    .. SEEALSO::

        :func:`~setup_intervals`
        :func:`~exists_vector`

    EXAMPLES::

        sage: from elementary_vectors import construct_vector, setup_intervals
        sage: M = matrix([1, 1, 0])
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    First, we consider closed intervals::

        sage: I = setup_intervals(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: construct_vector(M, I)
        (5, 5, 0)

    Next, we take open intervals. This time, there is no solution::

        sage: I = setup_intervals(lower_bounds, upper_bounds, False, False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
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
        sage: I = setup_intervals(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: construct_vector(M, I)
        (2, 5, 7, 5)
    
    TESTS:

    This example shows the case ``upper_bounds_closed == 1``::

        sage: M = matrix([[0, 0, 0, 0, 0, 0], [-1, -1, 0, -1, 1, 0], [-1, 1, 0, -2, -1, -2]])
        sage: I = setup_intervals(
        ....:     [-2, -1/2, -1, -1/4, -1/6, 0],
        ....:     [3/4, 1, 1, 0, 12, 1],
        ....:     [False, True, True, True, False, True],
        ....:     [False, True, False, False, False, True]
        ....: )
        sage: construct_vector(M, I)
        (-1/4, -1/4, 0, -1/4, 1/4, 0)
    
    zero matrix::
    
        sage: M = zero_matrix(QQ, 1, 5)
        sage: I = setup_intervals(
        ....:     [-1, -1, -1, -1, -1],
        ....:     [1, 1, 1, 1, 1],
        ....:     False,
        ....:     True
        ....: )
        sage: construct_vector(M, I)
        (0, 0, 0, 0, 0)
    """
    if not exists_vector(evs if evs else M, intervals):
        raise ValueError("There is no solution.")

    def rec(M, intervals):
        r"""Recursive call."""
        length = M.ncols()
        rank = M.rank()
        if rank == length:
            return simplest_vector_in_intervals(intervals)
        if rank == length - 1:
            # The kernel of ``M`` has one dimension.
            return construct_normal_vector(M.right_kernel_matrix().row(0), intervals)
        if rank == 1:
            for row in M.rows():
                if row != 0:
                    return multiple_in_intervals(row, intervals)
        if rank == 0:
            return zero_vector(M.base_ring(), M.ncols())

        def coefficient(k):
            r"""Construct ``x_k`` recursively"""
            # projection to k-th component
            M_bar = M.delete_columns([k])
            intervals_bar = intervals[:k] + intervals[k + 1:]
            xk_bar = rec(M_bar, intervals_bar)
            if hasattr(xk_bar, "simplify_full"):
                xk_bar = xk_bar.simplify_full()

            # solve linear system to project back
            try:
                x_k = M_bar.solve_left(xk_bar) * M # fails for certain matrices involving roots
            except ValueError:
                x_k = solve_left(M_bar, xk_bar) * M

            return x_k.simplify_full() if hasattr(x_k, "simplify_full") else x_k

        x = [coefficient(k) for k in range(rank + 2)]

        # use vectors in x to construct vector
        A = matrix([xk - x[0] for xk in x[1:]])
        weights = list(A.T.right_kernel_matrix().row(0))
        weights = [-sum(weights)] + weights
        v = sum(a * x_k for a, x_k in zip(weights, x) if a > 0) / sum(a for a in weights if a > 0)
        return v.simplify_full() if hasattr(v, "simplify_full") else v

    return rec(M, intervals)

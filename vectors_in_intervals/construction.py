"""Constructing vectors with components in intervals."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.calculus.var import var
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.infinity import Infinity
from sage.sets.real_set import RealSet
from sage.symbolic.relation import solve

from .existence import exists_orthogonal_vector, exists_vector
from .utility import simplest_element_in_interval, interval_from_bounds, solve_left
from elementary_vectors import elementary_vectors
from sign_vectors import sign_vector


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

        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.construction import *
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: simplest_vector_in_intervals(I)
        (2, 5, 0)
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, False, False)
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

        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.construction import *
        sage: I = intervals_from_bounds([0, 0, -1/4, -oo], [2, oo, 1/5, 1/9], False, [True, False, False, True])
        sage: v = vector([1, 5, -2, 0])
        sage: multiple_in_intervals_candidates(v, I)
        (0, 1/8)
        sage: v = vector([0, 5, -2, 0])
        sage: multiple_in_intervals_candidates(v, I)
        {}
        sage: I = intervals_from_bounds([-2, -1/4, -7], [0, 1, 0], [True, False, False], [False, False, True])
        sage: v = vector([-2, -2, -1])
        sage: multiple_in_intervals_candidates(v, I)
        (0, 1/8)
        sage: v = vector([0, -2, 1])
        sage: I = intervals_from_bounds([-2, 0, -1], [1, 1/3, 0], True, False)
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

    return interval_from_bounds(a_inf, a_sup, lower_bound_closed, upper_bound_closed)


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

        sage: from vectors_in_intervals import *
        sage: from vectors_in_intervals.construction import *
        sage: I = intervals_from_bounds([0, 0, -1/4, -oo], [2, oo, 1/5, 1/9], False, [True, False, False, True])
        sage: v = vector([1, 5, -2, 0])
        sage: multiple_in_intervals(v, I)
        (1/9, 5/9, -2/9, 0)
    """
    return simplest_element_in_interval(multiple_in_intervals_candidates(v, intervals)) * v


def vector_from_sign_vector(sv, data):
    r"""
    Find a vector in the row space of a matrix that has given signs.

    INPUT:

    - ``sv`` -- a sign vector of length ``n``
    
    - ``data`` -- either a real matrix with ``n`` columns or a list of
                elementary vectors of length ``n``
    
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
        sage: vector_from_sign_vector(zero_sign_vector(4), M)
        (0, 0, 0, 0)
        sage: vector_from_sign_vector(sign_vector("+-+0"), M)
        (2, -2, 2, 0)
        sage: vector_from_sign_vector(sign_vector("+0+0"), M)
        (1, 0, 2, 0)
        sage: vector_from_sign_vector(sign_vector("+-0+"), M)
        (1, -2, 0, 2)
        sage: vector_from_sign_vector(sign_vector("+-0+"), M)
        (1, -2, 0, 2)
        sage: evs = elementary_vectors(M.right_kernel_matrix())
        sage: vector_from_sign_vector(sign_vector("+-0+"), evs)
        (1, -2, 0, 2)
        sage: vector_from_sign_vector(sign_vector("+0-0"), M)
        Traceback (most recent call last):
        ...
        ValueError: There is no solution.
    """
    if isinstance(data, list):
        evs = data
    else:
        evs = elementary_vectors(data.right_kernel_matrix())

    def both_signs(evs):
        for v in evs:
            yield v
            yield -v

    result = sum(v for v in both_signs(evs) if sign_vector(v) <= sv)

    if isinstance(result, int): # empty sum is 0 (no evs or no solution)
        if isinstance(data, list):
            try:
                result = data[0].parent().zero_vector()
            except IndexError:
                result = zero_vector(sv.length())
        else:
            result = zero_vector(data.base_ring(), data.ncols())

    if sign_vector(result) != sv:
        raise ValueError("There is no solution.")
    return result


def construct_orthogonal_vector(v, intervals):
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

        :func:`vectors_in_intervals.existence.exists_orthogonal_vector`

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: I = intervals_from_bounds([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1, 1, 1])
        sage: construct_orthogonal_vector(v, I)
        (0, 1, -1)

    We define another vector. This time, there is no solution::

        sage: v = vector([1, 1, -1])
        sage: construct_orthogonal_vector(v, I)
        Traceback (most recent call last):
        ...
        ValueError: There is no solution.

    Next, we consider open intervals::

        sage: I = intervals_from_bounds([0, 1, -1], [1, 2, 1], False, [True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1, 0, 1])
        sage: construct_orthogonal_vector(v, I)
        (1/2, 2, -1/2)

    We can even consider unbounded intervals::

        sage: I = intervals_from_bounds([0, 1, -oo], [oo, 2, -2], False, [True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1, 1, -1])
        sage: construct_orthogonal_vector(v, I)
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
            continue
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
        if entry == 0:
            continue
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
            if product == 0:
                continue
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

        :func:`vectors_in_intervals.existence.exists_vector`

    EXAMPLES::

        sage: from vectors_in_intervals import *
        sage: M = matrix([1, 1, 0])
        sage: lower_bounds = [2, 5, -1]
        sage: upper_bounds = [5, 6, 1]

    First, we consider closed intervals::

        sage: I = intervals_from_bounds(lower_bounds, upper_bounds)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: construct_vector(M, I)
        (5, 5, 0)

    Next, we take open intervals. This time, there is no solution::

        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, False, False)
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
        sage: I = intervals_from_bounds(lower_bounds, upper_bounds, lower_bounds_closed, upper_bounds_closed)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: construct_vector(M, I)
        (2, 5, 7, 5)
    
    TESTS:

    This example shows the case ``upper_bounds_closed == 1``::

        sage: M = matrix([[0, 0, 0, 0, 0, 0], [-1, -1, 0, -1, 1, 0], [-1, 1, 0, -2, -1, -2]])
        sage: I = intervals_from_bounds(
        ....:     [-2, -1/2, -1, -1/4, -1/6, 0],
        ....:     [3/4, 1, 1, 0, 12, 1],
        ....:     [False, True, True, True, False, True],
        ....:     [False, True, False, False, False, True]
        ....: )
        sage: construct_vector(M, I)
        (-1/4, -1/4, 0, -1/4, 1/4, 0)
    
    zero matrix::
    
        sage: M = zero_matrix(QQ, 1, 5)
        sage: I = intervals_from_bounds(
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
            return construct_orthogonal_vector(M.right_kernel_matrix().row(0), intervals)
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
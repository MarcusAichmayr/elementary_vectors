r"""
Finding vectors in intervals.

With this module, we can check whether there is a vector in a subspace such that
the components lie in given intervals.

There is also an algorithmic approach to construct such a vector.
"""

#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .functions import elementary_vectors
from sage.sets.real_set import RealSet
from sage.rings.infinity import Infinity
from sage.calculus.var import var
from sage.symbolic.relation import solve
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix


def setup_intervals(L, R, l=True, r=True):
    r"""
    Construct a list of intervals from lists of bounds.

    INPUT:

    - ``L`` -- a list of real values and infinity of length ``n``

    - ``R`` -- a list of real values and infinity of length ``n``

    - ``l`` -- a boolean (default: ``True``) or a list of booleans of length ``n``

    - ``r`` -- a boolean (default: ``True``) or a list of booleans of length ``n``

    OUTPUT:

    A list of ``RealSet`` objects of length ``n``.

    - ``L`` and ``R`` are the left and right interval values.
      If ``L[i]`` > ``R[i]``, those elements will be exchanged.

    - ``l`` and ``r`` determine the intervals.

    - The left (or right) interval half of the ``i``-th interval is

        - closed if ``l[i]`` (or ``r[i]``) is ``True`` (default).

        - open if ``l[i]`` (or ``r[i]``) is ``False``.

    - If ``l`` (or ``r``) is a boolean, then all left (or right) interval halves
      are considered closed if ``True`` (default) and open if ``False``.

    EXAMPLES::

        sage: from elementary_vectors.vectors_in_intervals import setup_intervals
        sage: L = [2,5,-1]
        sage: R = [5,6,1]

    By default, the intervals are closed::

        sage: setup_intervals(L, R)
        [[2, 5], [5, 6], [-1, 1]]

    We obtain open intervals if both ``l`` and ``r`` are false::

        sage: setup_intervals(L, R, l=False, r=False)
        [(2, 5), (5, 6), (-1, 1)]

    Mixed intervals are also possible::

        sage: setup_intervals(L, R, l=False)
        [(2, 5], (5, 6], (-1, 1]]
        sage: setup_intervals(L, R, l=[False, True, True], r=[True, True, False])
        [(2, 5], [5, 6], [-1, 1)]

    We can also specify unbounded intervals.
    Note that bounds at infinity are always open::

        sage: L = [-oo, 2, -oo]
        sage: R = [-5, oo, oo]
        sage: setup_intervals(L, R)
        [(-oo, -5], [2, +oo), (-oo, +oo)]

    Finite and empty intervals are represented as usual::

        sage: L = [0, 2]
        sage: R = [0, 2]
        sage: setup_intervals(L, R, l=[True, False], r=True)
        [{0}, {}]

    If the lower bound is greater than the upper bound, those bounds will be exchanged::

        sage: L = [1, 4]
        sage: R = [2, 3]
        sage: setup_intervals(L, R, l=False)
        [(1, 2], (3, 4]]
    """
    n = len(L)
    if len(R) != n:
        raise ValueError('``R`` should be a list of length ' + str(n) + '.')

    if l is True:
        l = [True]*n
    elif l is False:
        l = [False]*n
    elif len(l) != n:
        raise ValueError('``l`` should be a list of length ' + str(n) + '.')

    if r is True:
        r = [True]*n
    elif r is False:
        r = [False]*n
    elif len(r) != n:
        raise ValueError('``r`` should be a list of length ' + str(n) + '.')

    intervals = []
    for Li, Ri, li, ri in zip(L, R, l, r):
        # if Li is the upper bound oo, we want to exchange the bounds first.
        if Ri < Li:
            Li, Ri = (Ri, Li)
        if Li == -Infinity:
            li = False
        if Ri == Infinity:
            ri = False

        if li and ri:
            I = RealSet.closed(Li, Ri)
        elif (not li) and (not ri):
            I = RealSet.open(Li, Ri)
        elif li and (not ri):
            I = RealSet.closed_open(Li, Ri)
        else:
            I = RealSet.open_closed(Li, Ri)
        intervals.append(I)
    return intervals


def exists_vector(data, intervals, kernel=False, certificate=False):
    r"""
    Return whether a vector exists in the vector space determined by a matrix such that the components lie in given intervals.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                  elementary vectors of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    - ``kernel`` -- a boolean (default: ``False``)

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:

    - If ``data`` is a matrix, then the elementary vectors of ``M`` are computed.

        - If ``kernel`` is false, considers the vector space generated by the
          rows of the matrix ``M``. (default)
          In this case, the elementary vectors will lie in the kernel of ``M``.

        - If ``kernel`` is true, considers the vector space generated by the
          kernel of the matrix ``M``.
          In this case, the elementary vectors will lie in the row space of ``M``.

    - If ``data`` is a list of elementary vectors, then those will be used.
      In this case, the argument ``kernel`` will be ignored.

    - If ``certificate`` is ``False`` (default), returns a boolean.

    - If ``certificate`` is ``True`` and the result is false, then a list
      ``[False, v]`` will be returned. Here, ``v`` is an elementary vector of ``M``
      that certifies that there exists no vector.

    ALGORITHM:

    The underlying algorithm is based on Minty's Lemma. (see [Min74]_)

    .. [Min74] Minty, G. J.:
       „A `from scratch` proof of a theorem of Rockafellar and Fulkerson“.
       In: Mathematical Programming 7 (1974), pp. 368-375.

    .. SEEALSO::

        :func:`~setup_intervals`

    EXAMPLES::

        sage: from elementary_vectors import exists_vector, setup_intervals
        sage: M = matrix([1, 1, 0])
        sage: L = [2, 5, -1] # lower halves of the intervals
        sage: R = [5, 6, 1] # upper halves of the intervals

    First, we consider closed intervals::

        sage: I = setup_intervals(L, R)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: exists_vector(M, I)
        True

    Next, we take open intervals::

        sage: I = setup_intervals(L, R, l=False, r=False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: exists_vector(M, I)
        False

    Mixed intervals are also possible::

        sage: l = [True, True, False]
        sage: r = [False, True, True]
        sage: I = setup_intervals(L, R, l, r)
        sage: I
        [[2, 5), [5, 6], (-1, 1]]
        sage: exists_vector(M, I)
        False

    Finally, we consider unbounded intervals::

        sage: M = matrix([[1, 0, 1, 0], [0, 1, 1, 1]])
        sage: L = [2, 5, 0, -oo]
        sage: R = [5, oo, 8, 5]
        sage: l = [True, True, False, False]
        sage: r = [False, False, False, True]
        sage: I = setup_intervals(L, R, l, r)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: exists_vector(M, I)
        True
    """
    if isinstance(data, list):
        evs = data
    else:
        evs = elementary_vectors(data, kernel=not kernel)

    for I in intervals:
        if I.is_empty():
            return False

    for v in evs:
        UB = 0
        LB = 0
        UBeq = True
        LBeq = True
        for vi, I in zip(v, intervals):
            # multiplication following 0 * oo = 0
            a, b = (0, 0) if vi == 0 else (vi*I.inf(), vi*I.sup())

            if a == b and vi != 0:
                LB += a
                UB += b
                LBeq &= (I.inf() in I) and (I.sup() in I)
                UBeq &= (I.inf() in I) and (I.sup() in I)
            elif a < b:
                LB += a
                UB += b
                LBeq &= I.inf() in I
                UBeq &= I.sup() in I
            elif a > b:
                LB += b
                UB += a
                LBeq &= I.sup() in I
                UBeq &= I.inf() in I

        if (LBeq and LB > 0) or (not LBeq and LB >= 0) or (UBeq and UB < 0) or (not UBeq and UB <= 0):
            return [False, v] if certificate else False

    return True


def construct_normal_vector(v, intervals):
    r"""
    Construct a vector that is normal on a given vector and lies in the specified intervals.

    INPUT:

    - ``v`` -- a vector of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    OUTPUT:

    Returns a vector ``z`` such that the scalar product of ``z`` and ``v`` is zero
    and each component of ``z`` lies in the respective interval of the list ``intervals``.
    If no such vector exists, raises a ``ValueError`` instead.

    .. SEEALSO::

        :func:`~setup_intervals`
    """
    if not exists_vector([v], intervals):
        raise ValueError("There is no solution.")
    z_min = []
    z_max = []
    var('eps')
    var('lam')
    open_intervals = False
    unbounded = False
    eps_values = []
    lam_values = []
    for vk, I in zip(v, intervals):
        if vk == 0:
            z_min.append(I.an_element())
            z_max.append(I.an_element())
        else:
            if not I.is_closed():
                if I.sup() != Infinity and I.inf() != -Infinity:
                    eps_values.append((I.sup() - I.inf())/2)
                open_intervals = True
            if I.inf() == -Infinity or I.sup() == Infinity:
                if I.inf() != -Infinity:
                    lam_values.append(abs(I.inf()))
                if I.sup() != Infinity:
                    lam_values.append(abs(I.sup()))
                unbounded = True
            l = (I.inf() + (0 if I.inf() in I else eps)) if I.inf() != -Infinity else -lam
            r = (I.sup() + (0 if I.sup() in I else -eps)) if I.sup() != Infinity else lam

            if vk > 0:
                z_min.append(l)
                z_max.append(r)
            else:
                z_min.append(r)
                z_max.append(l)
    z_min = vector(z_min)
    z_max = vector(z_max)

    if open_intervals:
        for z_m in [z_min, z_max]:
            try:
                if eps in (z_m*v).variables():
                    sol = solve(z_m*v, eps, solution_dict=True)[0][eps]
                    eps_values.append(sol)
            except AttributeError:
                pass
        eps_min = min(eps_values)
        # eps might not occur in both z_min and z_max
        try:
            if lam in eps_min.variables():
                lam_values.append(solve(eps_min, lam, solution_dict=True)[0][lam])
        except AttributeError:
            pass
        try:
            z_min = z_min(eps=eps_min)
        except TypeError:
            pass
        try:
            z_max = z_max(eps=eps_min)
        except TypeError:
            pass
    if unbounded:
        if z_min*v != 0:
            try:
                lam_values.append(solve(z_min*v, lam, solution_dict=True)[0][lam])
            except (IndexError, TypeError):
                pass
        if z_max*v != 0:
            try:
                lam_values.append(solve(z_max*v, lam, solution_dict=True)[0][lam])
            except (IndexError, TypeError):
                pass

        lam_max = max(lam_values) + 1
        try:
            z_min = z_min(lam=lam_max)
        except TypeError:
            pass
        try:
            z_max = z_max(lam=lam_max)
        except TypeError:
            pass

    a = v*z_min
    b = v*z_max
    if a == 0:
        return z_min
    if b == 0:
        return z_max

    z = (b*z_min - a*z_max)/(b - a) # convex combination lies in the intervals

    return z


def construct_vector(M, intervals):
    r"""
    Return a vector of a given vectorspace such that the components lie in given intervals.

    INPUT:

    - ``M`` -- a matrix with ``n`` columns

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    OUTPUT:

    Returns a vector in the rowspace of ``M`` such that each component lies
    in the respective interval of the list ``intervals``.
    If no such vector exists, raises a ``ValueError`` instead.

    .. SEEALSO::

        :func:`~setup_intervals`
    """
    if not exists_vector(M, intervals):
        raise ValueError("There is no solution.")

    def rec(M, intervals):
        r"""Recursive call."""
        n = M.ncols()
        r = M.rank()
        if r == n:
            return vector(I.an_element() for I in intervals)
        elif r == n - 1:
            # The kernel of ``M`` has one dimension. Hence, there is only one elementary vector.
            y = elementary_vectors(M, kernel=True)[0]
            return construct_normal_vector(y, intervals)
        else:
            # construct n vectors x1, ..., xn
            x = []
            for k in range(n):
                # projection to k-th component
                M_bar = M.delete_columns([k])
                intervals_bar = intervals[:k] + intervals[k+1:]
                xk_bar = rec(M_bar, intervals_bar)
                # solve linear system to project back
                a = M_bar.solve_left(xk_bar)
                xk = a*M
                x.append(xk)

            # use vectors in x to construct vector
            A = matrix([xk - x[0] for xk in x[1:]])
            a = list(A.T.right_kernel_matrix().row(0)) # dim == 1?
            a = [-sum(a)] + a
            sol = sum([ak*xk for ak, xk in zip(a, x) if ak > 0]) / sum(ak for ak in a if ak > 0)

            return sol

    return rec(M, intervals)

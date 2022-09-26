r"""
Finding vectors in intervals.

With this module, we can check whether there is a vector in a subspace such that
the components lie in given intervals.

There is also an algorithmic approach to construct such a vector.

EXAMPLES:

First, we load the functions from the package::

    sage: from elementary_vectors import *

The package offers the function :func:`~setup_intervals` that helps us creating lists of intervals.
This function takes two lists as input, the corresponding elements in those lists determine the intervals::

    sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
    sage: I
    [[0, 1], [1, 2], {-1}]

Next, we define a vector::

    sage: v = vector([1,1,1])

Is there a normal vector, such that the components lie in the intervals defined above?
We call :func:`~exists_normal_vector` to answer this question::

    sage: exists_normal_vector(v, I)
    True

This vector can be constructed by :func:`~construct_normal_vector`::

    sage: construct_normal_vector(v, I)
    (0, 1, -1)

We define another vector. This time, there is no solution::

    sage: v = vector([1,1,-1])
    sage: exists_normal_vector(v, I)
    False
    sage: construct_normal_vector(v, I)
    Traceback (most recent call last):
    ...
    ValueError: There is no solution.

Next, we consider open intervals::

    sage: I = setup_intervals([0, 1, -1], [1, 2, 1], l=False, r=[True, True, False])
    sage: I
    [(0, 1], (1, 2], (-1, 1)]
    sage: v = vector([1,0,1])
    sage: exists_normal_vector(v, I)
    True
    sage: construct_normal_vector(v, I)
    (1/2, 2, -1/2)

We can even consider unbounded intervals::

    sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], l=False, r=[True, True, False])
    sage: I
    [(0, +oo), (1, 2], (-oo, -2)]
    sage: v = vector([-1,1,-1])
    sage: exists_normal_vector(v, I)
    True
    sage: construct_normal_vector(v, I)
    (5, 2, -3)

The most important functions of this module are :func:`~exists_vector` and :func:`~construct_vector`.
Given a matrix ``M`` and a list of intervals, we want to examine whether there exists a vector in the rowspace of ``M``,
such that the components lie in the given intervals::

    sage: M = matrix([1, 1, 0])
    sage: L = [2, 5, -1]
    sage: R = [5, 6, 1]

First, we consider closed intervals::

    sage: I = setup_intervals(L, R)
    sage: I
    [[2, 5], [5, 6], [-1, 1]]
    sage: exists_vector(M, I)
    True
    sage: construct_vector(M, I)
    (5, 5, 0)

Next, we take open intervals. This time, there is no solution::

    sage: I = setup_intervals(L, R, l=False, r=False)
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
    sage: L = [2, 5, 0, -oo]
    sage: R = [5, oo, 8, 5]
    sage: l = [True, True, False, False]
    sage: r = [False, False, False, True]
    sage: I = setup_intervals(L, R, l, r)
    sage: I
    [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
    sage: exists_vector(M, I)
    True
    sage: construct_vector(M, I)
    (2, 5, 7, 5)
"""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
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
      If ``L[i] > R[i]``, those elements will be exchanged.

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

    return [setup_interval(Li, Ri, li, ri) for Li, Ri, li, ri in zip(L, R, l, r)]


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
        sage: L = [2,5,-1]
        sage: R = [5,6,1]
        sage: I = setup_intervals(L, R)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: simplest_vector_in_intervals(I)
        (2, 5, 0)
        sage: I = setup_intervals(L, R, l=False, r=False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: simplest_vector_in_intervals(I)
        (3, 11/2, 0)
    """
    return vector(simplest_element_in_interval(I) for I in intervals)


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
    lower_closed = False
    upper_closed = False

    for vk, Ik in zip(v, intervals):
        if vk == 0:
            if not 0 in Ik:
                return RealSet()
        else:
            val_inf = Ik.inf() / vk
            val_sup = Ik.sup() / vk
            if vk > 0:
                if val_inf > a_inf:
                    a_inf = val_inf
                    lower_closed = Ik.inf() in Ik
                elif val_inf == a_inf:
                    lower_closed = Ik.inf() in Ik and lower_closed
                if val_sup < a_sup:
                    a_sup = val_sup
                    upper_closed = Ik.sup() in Ik
                elif val_sup == a_sup:
                    upper_closed = Ik.sup() in Ik and upper_closed
            else: # vk < 0
                if val_sup > a_inf:
                    a_inf = val_sup
                    lower_closed = Ik.sup() in Ik
                elif val_sup == a_inf:
                    lower_closed = Ik.sup() in Ik and lower_closed
                if val_inf < a_sup:
                    a_sup = val_inf
                    upper_closed = Ik.inf() in Ik
                elif val_inf == a_sup:
                    upper_closed = Ik.inf() in Ik and upper_closed
            if a_inf > a_sup:
                return RealSet()

    return setup_interval(a_inf, a_sup, lower_closed, upper_closed)


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


def exists_normal_vector(v, intervals):
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

        sage: from elementary_vectors import exists_normal_vector, setup_intervals
        sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1,1,1])
        sage: exists_normal_vector(v, I)
        True
        sage: v = vector([1,1,-1])
        sage: exists_normal_vector(v, I)
        False

    Next, we consider open intervals::

        sage: I = setup_intervals([0, 1, -1], [1, 2, 1], l=False, r=[True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1,0,1])
        sage: exists_normal_vector(v, I)
        True

    We can even consider unbounded intervals::

        sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], l=False, r=[True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1,1,-1])
        sage: exists_normal_vector(v, I)
        True

    TESTS::

        sage: I = setup_intervals([-oo, 0], [oo, 0], l=False, r=False)
        sage: I
        [(-oo, +oo), {}]
        sage: v = vector([1,0])
        sage: exists_normal_vector(v, I)
        False
        sage: v = vector([1])
        sage: exists_normal_vector(v, I)
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

    for vi, I in zip(v, intervals):
        if I.is_empty():
            return False
        if vi != 0:
            if I.is_finite():  # I consists of one point
                lower_bound += vi*I.an_element()
                upper_bound += vi*I.an_element()
            elif vi > 0:
                lower_bound += vi*I.inf()
                upper_bound += vi*I.sup()
                lower_bound_zero &= I.inf() in I
                upper_bound_zero &= I.sup() in I
            else:  # vi < 0
                lower_bound += vi*I.sup()
                upper_bound += vi*I.inf()
                lower_bound_zero &= I.sup() in I
                upper_bound_zero &= I.inf() in I

    if lower_bound > 0:
        return False
    if upper_bound < 0:
        return False
    if lower_bound == 0 and not lower_bound_zero:
        return False
    if upper_bound == 0 and not upper_bound_zero:
        return False
    return True


def exists_vector(data, intervals):
    r"""
    Return whether a vector exists in a given vector space such that the components lie in the specified intervals.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                  elementary vectors of length ``n``

    - ``intervals`` -- a list of ``n`` intervals (``RealSet``)

    OUTPUT:

    Return whether there exists a vector in the vector space determined by data
    such that the components lie in the specified intervals.

    - If ``data`` is a matrix, then the elementary vectors in the kernel of ``M`` are computed.

    - If ``data`` is a list of elementary vectors, then those will be used.

    ALGORITHM:

    The underlying algorithm is based on Minty's Lemma. (see [Min74]_)

    .. [Min74] Minty, G. J.:
       „A `from scratch` proof of a theorem of Rockafellar and Fulkerson“.
       In: Mathematical Programming 7 (1974), pp. 368-375.

    .. SEEALSO::

        :func:`~setup_intervals`
        :func:`~exists_normal_vector`

    EXAMPLES::

        sage: from elementary_vectors import exists_vector, setup_intervals
        sage: M = matrix([1, 1, 0])
        sage: L = [2, 5, -1]
        sage: R = [5, 6, 1]

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

    TESTS::

        sage: I = setup_intervals([0, 0], [1, 0], l=False)
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
        evs = elementary_vectors(data, kernel=True, generator=True)
    else:
        for interval in intervals:
            if interval.is_empty():
                return False
        evs = data
    
    checked_supports = set()
    for v in evs:
        s = tuple(v.support())
        if not s in checked_supports:
            checked_supports.add(s)
            if not exists_normal_vector(v, intervals):
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
        :func:`~exists_normal_vector`

    EXAMPLES::

        sage: from elementary_vectors import construct_normal_vector, setup_intervals
        sage: I = setup_intervals([0, 1, -1], [1, 2, -1])
        sage: I
        [[0, 1], [1, 2], {-1}]
        sage: v = vector([1,1,1])
        sage: construct_normal_vector(v, I)
        (0, 1, -1)

    We define another vector. This time, there is no solution::

        sage: v = vector([1,1,-1])
        sage: construct_normal_vector(v, I)
        Traceback (most recent call last):
        ...
        ValueError: There is no solution.

    Next, we consider open intervals::

        sage: I = setup_intervals([0, 1, -1], [1, 2, 1], l=False, r=[True, True, False])
        sage: I
        [(0, 1], (1, 2], (-1, 1)]
        sage: v = vector([1,0,1])
        sage: construct_normal_vector(v, I)
        (1/2, 2, -1/2)

    We can even consider unbounded intervals::

        sage: I = setup_intervals([0, 1, -oo], [oo, 2, -2], l=False, r=[True, True, False])
        sage: I
        [(0, +oo), (1, 2], (-oo, -2)]
        sage: v = vector([-1,1,-1])
        sage: construct_normal_vector(v, I)
        (5, 2, -3)
    """
    if not exists_normal_vector(v, intervals):
        raise ValueError("There is no solution.")

    # construct z_min and z_max
    z_min = []
    z_max = []
    var('eps')
    var('lam')

    for vk, I in zip(v, intervals):
        if vk == 0:
            z_min.append(simplest_element_in_interval(I))
            z_max.append(simplest_element_in_interval(I))
        else:
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

    # find candidates for eps and lam
    eps_candidates = []
    lam_candidates = []
    for vk, I in zip(v, intervals):
        if vk != 0:
            if I.inf() == -Infinity and I.sup() == Infinity:  # (-oo, oo)
                lam_candidates.append(0)
            elif I.inf() == -Infinity:
                if I.sup() in I:  # (-oo, b]
                    lam_candidates.append(-I.sup())
                else:  # (-oo, b)
                    lam_candidates.append(1 - I.sup())
                    eps_candidates.append(1)
            elif I.sup() == Infinity:
                if I.inf() in I:  # [a, oo)
                    lam_candidates.append(I.inf())
                else:  # (a, oo)
                    lam_candidates.append(1 + I.inf())
                    eps_candidates.append(1)
            elif not I.is_closed():
                if I.sup() in I or I.inf() in I:  # [a, b) or (a, b]
                    eps_candidates.append(I.sup() - I.inf())
                else:  # (a, b)
                    eps_candidates.append((I.sup() - I.inf())/2)

    # determine eps
    if eps_candidates:
        for eq in [z_min*v, z_max*v]:
            try:
                if eps in eq.variables() and lam not in eq.variables():
                    eps_candidates.append(solve(eq, eps, solution_dict=True)[0][eps])
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

    # determine lam
    if lam_candidates:
        for eq in [z_min*v, z_max*v]:
            if eq != 0:
                try:
                    lam_candidates.append(solve(eq, lam, solution_dict=True)[0][lam])
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

    a = v*z_min
    b = v*z_max
    if a == 0:
        return z_min
    if b == 0:
        return z_max

    z = (b*z_min - a*z_max)/(b - a)  # convex combination lies in the intervals

    return z


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
        sage: L = [2, 5, -1]
        sage: R = [5, 6, 1]

    First, we consider closed intervals::

        sage: I = setup_intervals(L, R)
        sage: I
        [[2, 5], [5, 6], [-1, 1]]
        sage: construct_vector(M, I)
        (5, 5, 0)

    Next, we take open intervals. This time, there is no solution::

        sage: I = setup_intervals(L, R, l=False, r=False)
        sage: I
        [(2, 5), (5, 6), (-1, 1)]
        sage: construct_vector(M, I)
        Traceback (most recent call last):
        ...
        ValueError: There is no solution.

    Finally, we consider unbounded intervals::

        sage: M = matrix([[1, 0, 1, 0], [0, 1, 1, 1]])
        sage: L = [2, 5, 0, -oo]
        sage: R = [5, oo, 8, 5]
        sage: l = [True, True, False, False]
        sage: r = [False, False, False, True]
        sage: I = setup_intervals(L, R, l, r)
        sage: I
        [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
        sage: construct_vector(M, I)
        (2, 5, 7, 5)
    
    TESTS:

    This example calls the case ``r == 1``::

        sage: M = matrix([[0, 0, 0, 0, 0, 0], [-1, -1, 0, -1, 1, 0], [-1, 1, 0, -2, -1, -2]])
        sage: I = setup_intervals(
        ....:     [-2, -1/2, -1, -1/4, -1/6, 0],
        ....:     [3/4, 1, 1, 0, 12, 1],
        ....:     [False, True, True, True, False, True],
        ....:     [False, True, False, False, False, True]
        ....: )
        sage: construct_vector(M, I)
        (-1/4, -1/4, 0, -1/4, 1/4, 0)
    """
    if not exists_vector(evs if evs else M, intervals):
        raise ValueError("There is no solution.")

    def rec(M, intervals):
        r"""Recursive call."""
        n = M.ncols()
        r = M.rank()
        if r == n:
            return simplest_vector_in_intervals(intervals)
        elif r == 1:
            for row in M.rows():
                if row != 0:
                    return multiple_in_intervals(row, intervals)
        elif r == n - 1:
            # The kernel of ``M`` has one dimension.
            y = M.right_kernel_matrix().row(0)
            return construct_normal_vector(y, intervals)
        else:
            def comp_xk(k):
                r"""
                Recursively construct ``x_k``.
                """
                # projection to k-th component
                M_bar = M.delete_columns([k])
                intervals_bar = intervals[:k] + intervals[k+1:]
                xk_bar = rec(M_bar, intervals_bar)
                if hasattr(xk_bar, "simplify_full"):
                    xk_bar = xk_bar.simplify_full()
                # solve linear system to project back
                
                try:
                    a = M_bar.solve_left(xk_bar) # does not work for certain matrices involving roots
                except ValueError:
                    a = solve_left(M_bar, xk_bar)
                xk = a*M
                return xk.simplify_full() if hasattr(xk, "simplify_full") else xk

            x = [comp_xk(k) for k in range(r + 2)]

            # use vectors in x to construct vector
            A = matrix([xk - x[0] for xk in x[1:]])
            a = list(A.T.right_kernel_matrix().row(0))
            a = [-sum(a)] + a
            v = sum(ak*xk for ak, xk in zip(a, x) if ak > 0) / sum(ak for ak in a if ak > 0)
            return v.simplify_full() if hasattr(v, "simplify_full") else v

    return rec(M, intervals)

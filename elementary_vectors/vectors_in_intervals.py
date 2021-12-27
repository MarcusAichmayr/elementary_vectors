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


def setup_intervals(L, R, l=True, r=True):
    r"""
    Construct a list of intervals from lists of bounds.

    INPUT:

    - ``L`` -- a list of real values (``-oo`` and ``oo`` are accepted) of length ``n``

    - ``R`` -- a list of real values (``-oo`` and ``oo`` are accepted) of length ``n``

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

    Finite and empty intervals received the usual representation::

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

    - ``intervals`` -- a list of intervals (``RealSet``)

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

        - In this case, the argument ``kernel`` will be ignored.

    - ``L`` and ``R`` are the left and right interval values, respectively.

    - ``l`` and ``r`` determine the intervals.

    - The left (or right) interval half of the ``i``-th interval is

        - closed if ``l[i]`` (or ``r[i]``) is ``True`` (default).

        - open if ``l[i]`` (or ``r[i]``) is ``False``.

    - If ``l`` (or ``r``) is a boolean, then all left (or right) interval halves
      are considered closed if ``True`` (default) and open if ``False``.

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
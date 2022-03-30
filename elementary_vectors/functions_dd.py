r"""
Double description.

EXAMPLES:

We consider the following matrix::

    sage: M = matrix([[1,2,0,0,3],[0,1,-1,2,1]])
    sage: M
    [ 1  2  0  0  3]
    [ 0  1 -1  2  1]

Next, we compute the elementary vectors of this matrix::

    sage: from elementary_vectors import *
    sage: evs = elementary_vectors(M)
    sage: evs
    [(-2, 1, 1, 0, 0),
     (4, -2, 0, 1, 0),
     (-1, -1, 0, 0, 1),
     (0, 0, -2, -1, 0),
     (3, 0, -1, 0, -1),
     (-6, 0, 0, -1, 2),
     (0, 3, 1, 0, -2),
     (0, -6, 0, 1, 4)]

We compute the corresponding sign vectors::

    sage: from sign_vectors import *
    sage: svs = [sign_vector(v) for v in evs]
    sage: svs
    [(-++00), (+-0+0), (--00+), (00--0), (+0-0-), (-00-+), (0++0-), (0-0++)]

Now, we define a vector ``a``.
We will use this vector as a third row in our matrix::

    sage: a = vector([1, 1, 1, 0, 0])
    sage: Ma = M.insert_row(2, a)
    sage: Ma
    [ 1  2  0  0  3]
    [ 0  1 -1  2  1]
    [ 1  1  1  0  0]

That way, we describe a different vector space.
The corresponding elementary vectors can be computed as before::

    sage: evs_a = elementary_vectors(Ma)
    sage: evs_a
    [(-2, 1, 1, 0, 0), (-3, 3, 0, -1, -1), (-3, 0, 3, 1, 1), (0, -3, 3, 2, 2)]

Similarly, we obtain the following sign vectors::

    sage: [sign_vector(v) for v in evs_a]
    [(-++00), (-+0--), (-0+++), (0-+++)]

A different approach is to use the double description method.
First, we compute two lists of sign vectors::

    sage: from elementary_vectors.functions_dd import *
    sage: E0, Ep = dd_input(M, svs, a)
    sage: E0
    [(-++00)]
    sage: Ep
    [(+-0+0), (++00-), (00++0), (+0-0-), (+00+-), (0++0-), (0+0--)]

Then, we use the computed lists, to compute the new list of sign vectors
by applying double description::

    sage: dd(E0, Ep)
    [(-++00), (+-0++), (-0+++), (0-+++)]

There, is also a convenient command that computed this list of sign vectors::

    sage: double_description(M, a)
    [(-++00), (+-0++), (-0+++), (0-+++)]
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

from sign_vectors import sign_vector
from sign_vectors.utility import adjacent
from sage.combinat.combination import Combinations
from .reductions import reduce_vectors
from .utility import vector_from_matrix
from .functions import elementary_vectors
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix


def dd_input(M, svs, a):
    r"""
    INPUT:

    - ``M`` -- a matrix

    - ``svs`` -- a list of sign vectors

    - ``a`` -- a vector

    OUTPUT:
    a tuple ``(E0, Ep)``
    where ``E0`` are sign vectors such that the corresponding elementary vectors ``v`` satisfy ``a v = 0``
    and ``Ep`` are sign vectors such that the corresponding elementary vectors ``v`` satisfy ``a v > 0``.

    .. SEEALSO::

        :func:`~dd`
        :func:`~double_description`
        :func:`elementary_vectors.functions.elementary_vectors`
    """
    E0 = []
    Ep = []
    A = sign_vector(a)
    for X in svs:
        if X.disjoint_support(A):
            E0.append(X)
        elif X.is_harmonious(A):
            Ep.append(X)
        elif X.is_harmonious(-A):
            Ep.append(-X)
        else:
            x = vector_from_matrix(M, X.support())
            s = a*x
            if s == 0:
                E0.append(X)
            else:
                Ep.append(sign_vector(x if s > 0 else -x))

    return E0, Ep


def dd(E0, Ep, **kwargs):
    r"""
    Compute the sign vectors corresponding to given lists of sign vectors.

    INPUT:

    - ``E0`` -- a list of sign vectors

    - ``Ep`` -- a list of sign vectors

    OUTPUT:
    TODO

    .. SEEALSO::

        :func:`~dd_input`
        :func:`~double_description`
        :func:`~elementary_vectors`
    """
    out = []

    def check_pair(X, Y):
        if X.is_harmonious(Y):
            if not any(all((e in X.support() or e in Y.support()) for e in v.support()) for v in E0):
                if adjacent(X, Y, Ep) and adjacent(-X, -Y, Ep):
                    return True
        return False

    out = [X & -Y for X, Y in Combinations(Ep, 2) if check_pair(X, -Y)]

    return E0 + reduce_vectors(out, cancel_factors=False, **kwargs)


def double_description(M, a):
    r"""
    INPUT:

    - ``M`` -- a matrix

    - ``a`` -- a vector

    OUTPUT:
    a list of elementary vectors ``v`` of ``data`` satisfying ``a v = 0``.

    .. SEEALSO::

        :func:`~dd_input`
        :func:`~dd`
        :func:`elementary_vectors.functions.elementary_vectors`

    EXAMPLES:

    We consider the following matrix::

        sage: M = matrix([[1,2,0,0,3],[0,1,-1,2,1]])
        sage: M
        [ 1  2  0  0  3]
        [ 0  1 -1  2  1]

    Next, we compute the elementary vectors of this matrix::

        sage: from elementary_vectors import *
        sage: evs = elementary_vectors(M)
        sage: evs
        [(-2, 1, 1, 0, 0),
         (4, -2, 0, 1, 0),
         (-1, -1, 0, 0, 1),
         (0, 0, -2, -1, 0),
         (3, 0, -1, 0, -1),
         (-6, 0, 0, -1, 2),
         (0, 3, 1, 0, -2),
         (0, -6, 0, 1, 4)]

    We compute the corresponding sign vectors::

        sage: from sign_vectors import *
        sage: svs = [sign_vector(v) for v in evs]
        sage: svs
        [(-++00), (+-0+0), (--00+), (00--0), (+0-0-), (-00-+), (0++0-), (0-0++)]

    Now, we define a vector ``a``.
    We will use this vector as a third row in our matrix::

        sage: a = vector([1, 1, 1, 0, 0])
        sage: Ma = M.insert_row(2, a)
        sage: Ma
        [ 1  2  0  0  3]
        [ 0  1 -1  2  1]
        [ 1  1  1  0  0]

    That way, we describe a different vector space.
    The corresponding elementary vectors can be computed as before::

        sage: evs_a = elementary_vectors(Ma)
        sage: evs_a
        [(-2, 1, 1, 0, 0), (-3, 3, 0, -1, -1), (-3, 0, 3, 1, 1), (0, -3, 3, 2, 2)]

    Similarly, we obtain the following sign vectors::

        sage: [sign_vector(v) for v in evs_a]
        [(-++00), (-+0--), (-0+++), (0-+++)]

    A different approach is to use the double description method::

        sage: double_description(M, a)
        [(-++00), (+-0++), (-0+++), (0-+++)]

    The output is identical up to multiplies.
    """
    evs = elementary_vectors(M)
    svs = [sign_vector(v) for v in evs]
    E0, Ep = dd_input(M, svs, a)
    return dd(E0, Ep)


def cocircuits_iterative(ai):
    r"""
    Compute the sign vectors corresponding to given vectors.

    INPUT:

    - ``ai`` -- a list of vectors

    OUTPUT:
    Compute iteratively the sign vectors corresponding to the elementary vectors
    determined by the given vectors ``ai``.

    .. SEEALSO::

        :func:`~dd_input`
        :func:`~dd`
        :func:`~double_description`
        :func:`elementary_vectors.functions.elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors.functions_dd import cocircuits_iterative
        sage: ai = [vector([1,-2,0,2,2]), vector([0,1,4,4,1]), vector([1,0,-1,0,0])]
        sage: cocircuits_iterative(ai)
        [(0+0-+), (+-+-0), (+0+-+), (+-+0-)]
    """
    n = ai[0].length()
    ring = ai[0].base_ring()
    M = matrix(ring, 0, n)
    svs = [sign_vector(v) for v in identity_matrix(n)]
    for a in ai:
        E0, Ep = dd_input(M, svs, a)
        svs = dd(E0, Ep)
        M = matrix(list(M) + [a])
    return svs

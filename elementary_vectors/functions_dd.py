r"""
Double description

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
    sage: Ma = M.stack(a)
    sage: Ma
    [ 1  2  0  0  3]
    [ 0  1 -1  2  1]
    [ 1  1  1  0  0]

That way, we describe a different vector space.
The corresponding elementary vectors can be computed as before::

    sage: evs_a = elementary_vectors(Ma)
    sage: evs_a
    [(-4, 2, 2, 0, 0), (-6, 6, 0, -2, -2), (-6, 0, 6, 2, 2), (0, -6, 6, 4, 4)]

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
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.combination import Combinations
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sign_vectors import sign_vector
from sign_vectors.utility import adjacent

from .functions import elementary_vectors
from .reductions import reduce_vectors
from .utility import kernel_vector_support_given


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
    for X in svs:
        s = determine_sign(X, a, M)
        if s == 0:
            E0.append(X)
        else:
            Ep.append(X if s == 1 else -X)
    return E0, Ep


def determine_sign(X, a, M=None):
    r"""
    Determine the sign of the corresponding scalar product.

    INPUT:

    - ``X`` -- a sign vector

    - ``a`` -- a vector

    - ``M`` -- a matrix (default: None)

    OUTPUT:
    First, the vector ``v`` corresponding to the sign vector ``sv`` is computed.
    Then, the sign of the scalar product ``a v`` is returned.
    If the sign can not be determined, the matrix ``M`` is used to find an appropriate sign vector for ``X``.
    In this case, if ``M`` is not specified, an exception is raised.

    .. NOTE::

        This might not work if ``X`` is not a cocircuit of
        the oriented matroid corresponding to the kernel of ``M``.

    EXAMPLES::

        sage: from sign_vectors import *
        sage: from elementary_vectors.functions_dd import *
        sage: M = matrix([[1, 1, 2, 3], [2, -1, 0, 0]])
        sage: X = sign_vector("++-0")
        sage: a = vector([1, 0, 0, 0])
        sage: determine_sign(X, a, M)
        1
        sage: a = vector([0, 0, 0, 1])
        sage: determine_sign(X, a)
        0
        sage: a = vector([1, 1, 0, 0])
        sage: determine_sign(X, a)
        1
        sage: a = vector([1, -1, 0, 0])
        sage: determine_sign(X, a, M)
        -1
        sage: determine_sign(-X, a, M)
        1
    """
    if X.disjoint_support(a):
        return 0
    elif X.is_harmonious(a):
        return 1
    elif X.is_harmonious(-a):
        return -1
    else:
        if M is None:
            raise ValueError("Sign could not be determined. Pass a suitable matrix to determine the sign.")
        x = kernel_vector_support_given(M, X.support())
        if sign_vector(x) != X:
            x = -x
        s = a*x
        if s == 0:
            return 0
        else:
            return 1 if s > 0 else -1


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

    out = (X & -Y for X, Y in Combinations(Ep, 2) if check_pair(X, -Y))

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
        sage: Ma = M.stack(a)
        sage: Ma
        [ 1  2  0  0  3]
        [ 0  1 -1  2  1]
        [ 1  1  1  0  0]

    That way, we describe a different vector space.
    The corresponding elementary vectors can be computed as before::

        sage: evs_a = elementary_vectors(Ma)
        sage: evs_a
        [(-4, 2, 2, 0, 0), (-6, 6, 0, -2, -2), (-6, 0, 6, 2, 2), (0, -6, 6, 4, 4)]

    Similarly, we obtain the following sign vectors::

        sage: [sign_vector(v) for v in evs_a]
        [(-++00), (-+0--), (-0+++), (0-+++)]

    A different approach is to use the double description method::

        sage: double_description(M, a)
        [(-++00), (+-0++), (-0+++), (0-+++)]

    The output is identical up to multiplies.
    """
    evs = elementary_vectors(M)
    svs = (sign_vector(v) for v in evs)
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
    svs = (sign_vector(v) for v in identity_matrix(n))
    for a in ai:
        E0, Ep = dd_input(M, svs, a)
        svs = dd(E0, Ep)
        M = matrix(list(M) + [a])
    return svs

r"""Computing elementary vectors."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import warnings
from sage.combinat.combination import Combinations
from sage.modules.free_module_element import vector, zero_vector
from sage.structure.element import get_coercion_model
from sage.arith.misc import gcd
from sign_vectors import sign_vector
from sign_vectors.utility import adjacent
from sage.symbolic.ring import SR
from .reductions import reduce_vectors, reduce_vector_using_equalities
from .utility import sign_determined
from sage.symbolic.assumptions import assume, forget


def elementary_vectors(data, dim=None, kernel=True, return_minors=False, ring=None, generator=False, **kwargs):
    r"""
    Compute elementary vectors of a subspace determined by the matrix or from a list of maximal minors.

    INPUT:

    - ``data`` -- a matrix or a list of maximal minors

    - ``dim`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.

    - ``kernel`` -- a boolean (default: ``True``)

    - ``return_minors`` -- a boolean (default: ``False``)

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``data``.

    - ``generator`` -- a boolean (default: ``False``)

    .. NOTE::

        Keyword arguments may be specified to apply certain reductions to the output.
        By default, all those reductions (like canceling factors) are applied.
        Possible keyword arguments are the same as in the function
        :func:`elementary_vectors.reductions.reduce_vectors`.

    OUTPUT:

    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``data``. (default)
      Otherwise, returns a list of elementary vectors lying in
      the row space of ``data``.
      This argument is ignored if ``data`` is not a matrix.

    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.
    
    - If ``generator`` is true, a generator of elementary vectors is returned.

    EXAMPLES::

        sage: from elementary_vectors import elementary_vectors
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors(M)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, kernel=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, kernel=False)
        [(0, 1, 2, 0, 0), (0, 1, 0, 2, 0), (2, 0, 0, 0, 2), (0, 0, 1, -1, 0)]
        sage: elementary_vectors(M, return_minors=True)
        [[(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)], [1, -1, 0, -2, 0, 0, 0, 1, -1, -2]]

    By default, the output is reduced, that is, each returned elementary vector
    has distinct support and common factors are canceled::

        sage: elementary_vectors(M)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, cancel_factors=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2)]
        sage: elementary_vectors(M, cancel_factors=False, reduce_support=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]

    Rational matrices are supported as well::

        sage: A = matrix([[1/2, 0, 0, 2/3], [0, 0, 5/2, -1]])
        sage: A
        [1/2   0   0 2/3]
        [  0   0 5/2  -1]
        sage: elementary_vectors(A)
        [(0, -15, 0, 0), (-20, 0, 6, 15)]
        sage: elementary_vectors(A, cancel_factors=False)
        [(0, -5/4, 0, 0), (-5/3, 0, 1/2, 5/4)]

    We can also compute elementary vectors over a finite field::

        sage: B = matrix(GF(7), [[1, 2, 3, 4, 0], [0, 5, 2, 3, 3]])
        sage: B
        [1 2 3 4 0]
        [0 5 2 3 3]
        sage: elementary_vectors(B)
        [(3, 5, 5, 0, 0),
         (0, 4, 0, 5, 0),
         (6, 4, 0, 0, 5),
         (1, 0, 4, 2, 0),
         (2, 0, 4, 0, 2),
         (5, 0, 0, 4, 3),
         (0, 2, 1, 0, 3),
         (0, 0, 5, 5, 1)]

    Variables are also supported::

        sage: var('c')
        c
        sage: C = matrix([[c, 0, 0, 2, 1, 1], [0, 1, 0, 1, 2, 1], [0, 0, 1, 1, 1, 2]])
        sage: C
        [c 0 0 2 1 1]
        [0 1 0 1 2 1]
        [0 0 1 1 1 2]
        sage: elementary_vectors(C)
        [(2, c, c, -c, 0, 0),
         (1, 2*c, c, 0, -c, 0),
         (1, c, 2*c, 0, 0, -c),
         (-1, c, 0, c, -c, 0),
         (-3, -c, 0, 2*c, 0, -c),
         (-1, -3*c, 0, 0, 2*c, -c),
         (3, 0, c, -2*c, c, 0),
         (1, 0, -c, -c, 0, c),
         (-1, 0, -3*c, 0, -c, 2*c),
         (4, 0, 0, -3*c, c, c),
         (0, 3, 1, 1, -2, 0),
         (0, 1, 3, 1, 0, -2),
         (0, -1, 1, 0, 1, -1),
         (0, 4, 0, 1, -3, 1),
         (0, 0, 4, 1, 1, -3)]

    Here, ``data`` is a list of maximal minors::

        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: m = M.minors(3); m
        [2, -2, 0, -4, 0, 0, 0, 2, -2, -4]
        sage: elementary_vectors(m, [3,5])
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors(m, M.dimensions())
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors(m, M.dimensions(), ring=QQ)
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]

    TESTS::

        sage: M = random_matrix(QQ, 0, 4)
        sage: elementary_vectors(M)
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    """
    if isinstance(data, list):
        if dim is None:
            raise TypeError("When computing elementary vectors from a list of " +
                            "maximal minors, the dimensions of the corresponding " +
                            "matrix are needed.")
        evs = elementary_vectors_from_minors(data, dim=dim, ring=ring, generator=generator, **kwargs)
        if return_minors:
            return [evs, data]
        else:
            return evs
    else:
        return elementary_vectors_from_matrix(data, kernel=kernel, return_minors=return_minors, ring=ring, generator=generator, **kwargs)


def elementary_vectors_from_matrix(M, kernel=True, return_minors=False, ring=None, generator=False, **kwargs):
    r"""
    Compute elementary vectors of a subspace determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    - ``return_minors`` -- a boolean (default: ``False``)

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``M``.

    .. NOTE::

        Keyword arguments may be specified to apply certain reductions to the output.
        By default, all those reductions (like canceling factors) are applied.
        Possible keyword arguments are the same as in the function
        :func:`elementary_vectors.reductions.reduce_vectors`.

    OUTPUT:

    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``M``. (default)

    - If ``kernel`` is false, returns a list of elementary vectors lying in
      the row space of ``M``.

    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.

    .. SEEALSO::

        :func:`~elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors.functions import elementary_vectors_from_matrix
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors_from_matrix(M)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_matrix(M, kernel=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_matrix(M, kernel=False)
        [(0, 1, 2, 0, 0), (0, 1, 0, 2, 0), (2, 0, 0, 0, 2), (0, 0, 1, -1, 0)]
        sage: elementary_vectors_from_matrix(M, return_minors=True)
        [[(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)], [1, -1, 0, -2, 0, 0, 0, 1, -1, -2]]
    """
    if kernel:
        try:
            ind = M.pivot_rows()  # would not work for polynomial matrices
            M = M.matrix_from_rows(ind)
            full_rank = True
        except (ArithmeticError, NotImplementedError):
            full_rank = False  # The matrix might not have full rank. In this case, we will obtain ``[]``.
    else:
        M = M.right_kernel_matrix()  # not implemented for polynomial matrices over ZZ (QQ does work.)

    if generator:
        return elementary_vectors_generator(M, **kwargs)

    n = M.ncols()
    r = M.nrows()
    m = M.minors(r)
    try:
        g = gcd(m)
        if g == 0:
            if not full_rank:
                warnings.warn('Could not determine rank of matrix. Result might be wrong.')
                return []
        try:
            if kwargs['cancel_factors']:
                m = [mi/g for mi in m]
        except KeyError:
            m = [mi/g for mi in m]
    except ValueError:
        pass
    if ring is None:
        ring = M.base_ring()
    evs = elementary_vectors_from_minors(m, [r, n], ring=ring, **kwargs)
    if return_minors:
        return [evs, m]
    else:
        return evs


def elementary_vectors_from_minors(m, dim, ring=None, generator=False, **kwargs):
    r"""
    Compute elementary vectors determined by given maximal minors of a matrix.

    INPUT:

    - ``m`` -- a list of maximal minors of a matrix

    - ``dim`` -- a tuple of the dimensions of the matrix corresponding to ``m``

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``m``.

    - ``generator`` -- a boolean (default: ``False``)

    .. NOTE::

        Keyword arguments may be specified to apply certain reductions to the output.
        By default, multiples and zero vectors are removed from the output.
        Factors are not canceled by default since this operation is less efficient.
        Possible keyword arguments are the same as in the function
        :func:`elementary_vectors.reductions.reduce_vectors`.

    OUTPUT:

    Uses the maximal minors ``m`` to compute the elementary vectors of the
    corresponding matrix.
    If ``generator`` is true, a generator of elementary vectors will be returned instead of a list.
    In this case, no reductions are applied.

    .. SEEALSO::

        :func:`~elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors.functions import elementary_vectors_from_minors
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: m = M.minors(3); m
        [2, -2, 0, -4, 0, 0, 0, 2, -2, -4]
        sage: elementary_vectors_from_minors(m, [3, 5])
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors_from_minors(m, M.dimensions())
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors_from_minors(m, M.dimensions(), ring=QQ)
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
    """
    r, n = dim
    dets = Combinations(n, r)
    if ring is None:
        ring = get_coercion_model().common_parent(*m)

    def ev_from_minors(I):
        r"""
        Compute the elementary vector corresponding to ``I``.

        INPUT:

        - ``I`` -- a list of indices
        """
        v = zero_vector(ring, n)
        for pos, k in enumerate(I):
            v[k] = (-1)**pos * m[dets.rank(i for i in I if i != k)]
        return v

    evs = (ev_from_minors(I) for I in Combinations(n, r + 1))
    if generator:
        return evs
    else:
        return reduce_vectors(evs, **kwargs)


def elementary_vectors_generator(M, **kwargs):
    r"""
    Compute a generator of elementary vectors.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    OUTPUT:
    
    a generator of all elementary vectors in the right kernel of ``M``.
    
    .. SEEALSO::

        :func:`~elementary_vectors`

    EXAMPLES::
    
        sage: from elementary_vectors.functions import elementary_vectors_generator
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: evs = elementary_vectors_generator(M)
        sage: next(evs)
        (0, 2, -1, -1, 0)
        sage: next(evs)
        (1, 0, 0, 0, -1)
    """
    m, n = M.dimensions()
    
    checked_supports = set()
    for I in Combinations(n, m + 1):
        M_I = M.matrix_from_columns(I)
        for v in elementary_vectors(M_I, **kwargs):
            l = list(v)
            for k in range(n):
                if not k in I:
                    l.insert(k, 0)
            v = vector(l)
            s = frozenset(v.support())
            if s not in checked_supports:
                checked_supports.add(s)
                yield(v)


def non_negative_vectors(L):
    r"""
    Return non-negative vectors.

    INPUT:

    - ``L`` -- a list of vectors

    OUTPUT:

    Return all vectors of ``L`` that are
    - non_negative in each component; or
    - negative in each component. Those will be multiplied by ``-1``; or
    - containing variables such that no opposing signs occur.

    EXAMPLES::

        sage: from elementary_vectors import non_negative_vectors
        sage: l = [vector([1, 1, 0, -1]), vector([0, 0, 0, 0]), vector([1, 0, 0, 1])]
        sage: l
        [(1, 1, 0, -1), (0, 0, 0, 0), (1, 0, 0, 1)]
        sage: non_negative_vectors(l)
        [(0, 0, 0, 0), (1, 0, 0, 1)]

    Now, we consider an example with a variable::

        sage: from elementary_vectors import elementary_vectors
        sage: var('a')
        a
        sage: A = matrix([[a, 0, 0, 0, 1], [0, 1, 0, 0, 1]])
        sage: evs = elementary_vectors(A, cancel_factors=True)
        sage: evs
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (-1, -a, 0, 0, a)]
        sage: non_negative_vectors(evs)
        ...
        UserWarning: Cannot determine sign of symbolic expression, returning 0 instead.
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (1, a, 0, 0, -a)]
        sage: assume(a > 0)
        sage: non_negative_vectors(evs)
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]
    """
    out = []
    for v in L:
        if sign_vector(v) >= 0:  # Use ``>=`` instead of ``>``, (0,0,x) -> (000) should be returned
            out.append(v)
        elif sign_vector(v) < 0:
            out.append(-v)
    return out


# TODO: should assumptions be an optional argument?
def positive_elementary_vectors(data, dim=None, kernel=True, return_minors=False, ring=None, **kwargs):
    r"""
    Compute positive elementary vectors.

    INPUT:

    - ``data`` -- a matrix or a list of maximal minors

    - ``dim`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.

    - ``kernel`` -- a boolean (default: ``True``)

    - ``return_minors`` -- a boolean (default: ``False``)

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``data``.

    .. NOTE::

        Keyword arguments may be specified to apply certain reductions to the output.
        By default, all those reductions (like canceling factors) are applied.
        Possible keyword arguments are the same as in the function
        :func:`elementary_vectors.reductions.reduce_vectors`.

    OUTPUT:

    The output is a list of pairs.
    Each pair consists of a list of assumptions and the corresponding positive elementary vectors of ``data``.

    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``data``. (default)
      Otherwise, returns a list of elementary vectors lying in
      the row space of ``data``.
      This argument is ignored if ``data`` is not a matrix.

    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.

    EXAMPLES::

        sage: from elementary_vectors import positive_elementary_vectors
        sage: A = matrix([[1, -1, 0]])
        sage: positive_elementary_vectors(A)
        [[[], [(1, 1, 0), (0, 0, 1)]]]
        sage: positive_elementary_vectors(A, return_minors=True)
        [[[], [(1, 1, 0), (0, 0, 1)], [1, -1, 0]]]
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]])
        sage: positive_elementary_vectors(M)
        [[[], []]]

        TODO: Do more examples also symbolic examples
    """
    kwargs["cancel_factors"] = False

    def rec(i, l, eq):
        if i < len(m):
            if not sign_determined(m[i]):
                a = SR(m[i])
                try:
                    expr = a > 0
                    assume(expr)
                    rec(i+1, l + [expr], eq)
                    forget(expr)
                except ValueError:
                    pass
                try:
                    expr = a < 0
                    assume(expr)
                    rec(i+1, l + [expr], eq)
                    forget(expr)
                except ValueError:
                    pass
                try:
                    expr = a == 0
                    assume(expr)
                    rec(i+1, l + [expr], eq + [expr, -a == 0])  # a.substitute(-a == 0) results in a instead of 0.
                    forget(expr)
                except ValueError:
                    pass
            else:
                rec(i+1, l, eq)
        else:
            m_eq = reduce_vector_using_equalities(vector(m), eq=eq)
            if m_eq == 0:  # minors are all zero
                warnings.warn('For ' + str(l) + ' all maximal minors are zero. There might be missing positive elementary vectors.')
#                 if isinstance(data, list):
#                     out.append([l, 'maximal minors are zero, could not compute elementary vectors for this branch'])
#                     return
#                 else:
#                     print('Todo')
#                     # substitute in data
#                     # compute positive_elementary_vectors from resulting matrix, eventuell auch l übergeben.
#                     # append each pair of result to out, also append assumptions l to first arguments
#                     pass
            # Do not overwrite ``evs`` here! It might be used in another branch.

            L = reduce_vectors(evs, eq=eq, **kwargs)  # TODO: this might result in an empty list. Instead, compute elementary vectors again if data is matrix
            if return_minors:
                out.append([l, non_negative_vectors(L), list(m_eq)])
            else:
                out.append([l, non_negative_vectors(L)])
            return

    evs, m = elementary_vectors(data, dim=dim, kernel=kernel, return_minors=True, ring=ring, **kwargs)
    if evs == []:
        return []
    out = []
    rec(0, [], [])
    return out

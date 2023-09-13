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
from sage.modules.free_module_element import zero_vector
from sage.structure.element import get_coercion_model

from sign_vectors import sign_vector

from .reductions import remove_multiples_generator


def elementary_vectors(data, dimensions=None, remove_multiples=True, ring=None, generator=False):
    r"""
    Compute elementary vectors of a subspace determined by a matrix or a list of maximal minors.

    INPUT:

    - ``data`` -- a matrix or a list of maximal minors

    - ``dimensions`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.

    - ``remove_multiples`` -- a boolean (default: ``False``)

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``data``.

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:

    If ``data`` is a matrix,
    returns a list of elementary vectors lying in the kernel of this matrix.

    - If ``generator`` is true, a generator of elementary vectors is returned.

    EXAMPLES::

        sage: from elementary_vectors import elementary_vectors
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors(M)
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]

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
        if dimensions is None:
            raise TypeError("When computing elementary vectors from a list of " +
                            "maximal minors, the dimensions of the corresponding " +
                            "matrix are needed.")
        return elementary_vectors_from_minors(data, dimensions=dimensions, ring=ring, remove_multiples=remove_multiples, generator=generator)
    if generator:
        return elementary_vectors_from_matrix(data, remove_multiples=remove_multiples, generator=True)
    return elementary_vectors_from_minors(data.minors(data.nrows()), data.dimensions(), ring=data.base_ring(), remove_multiples=remove_multiples)


def elementary_vectors_from_matrix(M, remove_multiples=True, generator=False):
    r"""
    Compute elementary vectors of a subspace determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``remove_multiples`` -- a boolean (default: ``False``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:

    Returns a list of elementary vectors lying in the kernel of ``M``.

    - If ``remove_multiples`` is true, the output will be reduced such that
      all elementary vectors have distinct support.
    
    - If ``generator`` is true, the output will be a generator object instead of a list.

    .. SEEALSO::

        :func:`~elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors.functions import elementary_vectors_from_matrix
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors_from_matrix(M)
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors_from_matrix(M, remove_multiples=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]
    """
    try:
        ind = M.pivot_rows()  # does not work for polynomial matrices
        M = M.matrix_from_rows(ind)
    except (ArithmeticError, NotImplementedError):
        warnings.warn('Could not determine rank of matrix. Result might be wrong.')

    rank, length = M.dimensions()
    minors = {}
    ring = M.base_ring()

    def ev_from_support(indices):
        r"""
        Return the elementary vector corresponding to ``indices``.

        INPUT:

        - ``indices`` -- a list of indices
        """
        element = zero_vector(ring, length)
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            try:
                minor = minors[indices_minor]
            except KeyError:
                minors[indices_minor] = M.matrix_from_columns(indices_minor).det()
                minor = minors[indices_minor]
            element[k] = (-1)**pos * minor
        return element

    evs = (ev_from_support(indices) for indices in Combinations(length, rank + 1))
    evs = (v for v in evs if v)
    if remove_multiples:
        evs = remove_multiples_generator(evs)

    if generator:
        return evs
    return list(evs)
#    return evs if generator else list(evs)


def elementary_vectors_from_minors(minors, dimensions, ring=None, remove_multiples=True, generator=False):
    r"""
    Compute elementary vectors determined by given maximal minors of a matrix.

    INPUT:

    - ``minors`` -- a list of maximal minors of a matrix

    - ``dim`` -- a tuple of the dimensions of the matrix corresponding to ``minors``

    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``minors``.

    - ``generator`` -- a boolean (default: ``False``)

    .. NOTE::

        Keyword arguments may be specified to apply certain reductions to the output.
        By default, multiples and zero vectors are removed from the output.
        Factors are not canceled by default since this operation is less efficient.
        Possible keyword arguments are the same as in the function
        :func:`elementary_vectors.reductions.reduce_vectors`.

    OUTPUT:

    Uses the maximal minors ``minors`` to compute the elementary vectors of the
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
    rank, length = dimensions
    if ring is None:
        ring = get_coercion_model().common_parent(*minors)

    minors_dict = {}
    for k, indices in enumerate(Combinations(length, rank)):
        minors_dict[tuple(indices)] = minors[k]

    def ev_from_support(indices):
        r"""
        Return the elementary vector corresponding to ``indices``.

        INPUT:

        - ``indices`` -- a list of indices
        """
        element = zero_vector(ring, length)
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            minor = minors_dict[indices_minor]
            element[k] = (-1)**pos * minor
        return element

    evs = (ev_from_support(indices) for indices in Combinations(length, rank + 1))
    evs = (v for v in evs if v)
    if remove_multiples:
        evs = remove_multiples_generator(evs)

    return evs if generator else list(evs)


def non_negative_vectors(vectors):
    r"""
    Return non-negative vectors.

    INPUT:

    - ``vectors`` -- a list of vectors

    OUTPUT:

    Return all vectors of ``vectors`` that are
    - non_negative in each component; or
    - negative in each component. Those will be multiplied by ``-1``; or
    - containing variables such that no opposing signs occur.

    EXAMPLES::

        sage: from elementary_vectors.functions import non_negative_vectors
        sage: l = [vector([1, 1, 0, -1]), vector([0, 0, 0, 0]), vector([1, 0, 0, 1])]
        sage: l
        [(1, 1, 0, -1), (0, 0, 0, 0), (1, 0, 0, 1)]
        sage: non_negative_vectors(l)
        [(0, 0, 0, 0), (1, 0, 0, 1)]

    Now, we consider an example with a variable::

        sage: from elementary_vectors import elementary_vectors
        sage: from elementary_vectors.reductions import *
        sage: var('a')
        a
        sage: A = matrix([[a, 0, 0, 0, 1], [0, 1, 0, 0, 1]])
        sage: evs = reduce_vectors(elementary_vectors(A), cancel_factors=True)
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
    result = []
    for element in vectors:
        if sign_vector(element) >= 0:  # ``>=`` instead of ``>``, (0,0,x) -> (000) should be returned
            result.append(element)
        elif sign_vector(element) < 0:
            result.append(-element)
    return result


# TODO: should assumptions be an optional argument?
# def positive_elementary_vectors(data, dim=None, kernel=True, return_minors=False, ring=None, **kwargs):
#     r"""
#     Compute positive elementary vectors.

#     INPUT:

#     - ``data`` -- a matrix or a list of maximal minors

#     - ``dim`` -- Not needed if ``data`` is a matrix.
#                  Otherwise, this is a list of dimensions of the matrix
#                  corresponding to the list of maximal minors ``data``.

#     - ``kernel`` -- a boolean (default: ``True``)

#     - ``return_minors`` -- a boolean (default: ``False``)

#     - ``ring`` -- Parent of the entries of the elementary vectors.
#                   By default, determine this from ``data``.

#     .. NOTE::

#         Keyword arguments may be specified to apply certain reductions to the output.
#         By default, all those reductions (like canceling factors) are applied.
#         Possible keyword arguments are the same as in the function
#         :func:`elementary_vectors.reductions.reduce_vectors`.

#     OUTPUT:

#     The output is a list of pairs.
#     Each pair consists of a list of assumptions and the corresponding positive elementary vectors of ``data``.

#     - If ``kernel`` is true, returns a list of elementary vectors lying in
#       the kernel of ``data``. (default)
#       Otherwise, returns a list of elementary vectors lying in
#       the row space of ``data``.
#       This argument is ignored if ``data`` is not a matrix.

#     - If ``return_minors`` is true, a list is returned where the first
#       element is the list of elementary vectors and the second element is
#       the list of maximal minors used to compute the elementary vectors.

#     EXAMPLES::

#         sage: from elementary_vectors import positive_elementary_vectors
#         sage: A = matrix([[1, -1, 0]])
#         sage: positive_elementary_vectors(A)
#         [[[], [(1, 1, 0), (0, 0, 1)]]]
#         sage: positive_elementary_vectors(A, return_minors=True)
#         [[[], [(1, 1, 0), (0, 0, 1)], [1, -1, 0]]]
#         sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]])
#         sage: positive_elementary_vectors(M)
#         [[[], []]]

#         TODO: Do more examples also symbolic examples
#     """
#     kwargs["cancel_factors"] = False

#     def rec(i, l, eq):
#         if i < len(m):
#             if not sign_determined(m[i]):
#                 a = SR(m[i])
#                 try:
#                     expr = a > 0
#                     assume(expr)
#                     rec(i+1, l + [expr], eq)
#                     forget(expr)
#                 except ValueError:
#                     pass
#                 try:
#                     expr = a < 0
#                     assume(expr)
#                     rec(i+1, l + [expr], eq)
#                     forget(expr)
#                 except ValueError:
#                     pass
#                 try:
#                     expr = a == 0
#                     assume(expr)
#                     rec(i+1, l + [expr], eq + [expr, -a == 0])  # a.substitute(-a == 0) results in a instead of 0.
#                     forget(expr)
#                 except ValueError:
#                     pass
#             else:
#                 rec(i+1, l, eq)
#         else:
#             m_eq = reduce_vector_using_equalities(vector(m), eq=eq)
#             if m_eq == 0:  # minors are all zero
#                 warnings.warn('For ' + str(l) + ' all maximal minors are zero. There might be missing positive elementary vectors.')
# #                 if isinstance(data, list):
# #                     result.append([l, 'maximal minors are zero, could not compute elementary vectors for this branch'])
# #                     return
# #                 else:
# #                     print('Todo')
# #                     # substitute in data
# #                     # compute positive_elementary_vectors from resulting matrix, eventuell auch l übergeben.
# #                     # append each pair of result to result, also append assumptions l to first arguments
# #                     pass
#             # Do not overwrite ``evs`` here! It might be used in another branch.

#             L = reduce_vectors(evs, eq=eq, **kwargs)  # TODO: this might result in an empty list. Instead, compute elementary vectors again if data is matrix
#             if return_minors:
#                 result.append([l, non_negative_vectors(L), list(m_eq)])
#             else:
#                 result.append([l, non_negative_vectors(L)])
#             return

#     evs, m = elementary_vectors(data, dim=dim, kernel=kernel, return_minors=True, ring=ring, **kwargs)
#     if evs == []:
#         return []
#     result = []
#     rec(0, [], [])
#     return result

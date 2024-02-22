r"""Computing elementary vectors."""

#############################################################################
#  Copyright (C) 2024                                                       #
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

from .reductions import reduce_vectors_support


def elementary_vectors(data, dimensions=None, remove_multiples=True, generator=False):
    r"""
    Compute elementary vectors of a subspace determined by a matrix or a list of maximal minors.

    INPUT:

    - ``data`` -- a matrix or a list of maximal minors

    - ``dimensions`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.

    - ``remove_multiples`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix given by ``data``.
    The element ``data`` can be a matrix or a list of maximal minors.
    In the latter case, dimensions are also required.

    - If ``remove_multiples`` is true, the output will be reduced such that
      all elementary vectors have distinct support.
    
    - If ``generator`` is true, the output will be a generator object instead of a list.

    .. SEEALSO::

        :func:`~elementary_vectors_from_matrix`
        :func:`~elementary_vectors_from_minors`

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

    Next, we determine the maximal minors of some matrix and use them as input.
    This is useful if we want to reuse the maximal minors or if they are given::

        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: m = M.minors(3); m
        [2, -2, 0, -4, 0, 0, 0, 2, -2, -4]
        sage: elementary_vectors(m, [3, 5])
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
        sage: elementary_vectors(m, M.dimensions())
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]

    TESTS::

        sage: elementary_vectors(random_matrix(QQ, 0, 4))
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        sage: elementary_vectors(matrix([[1, 0, 0], [0, 0, 0]]))
        [(0, -1, 0), (0, 0, -1)]
    """
    if isinstance(data, list):
        if dimensions is None:
            raise TypeError("When computing elementary vectors from a list of " +
                            "maximal minors, the dimensions of the corresponding " +
                            "matrix are needed.")
        return elementary_vectors_from_minors(data, dimensions=dimensions, remove_multiples=remove_multiples, generator=generator)
    return elementary_vectors_from_matrix(data, remove_multiples=remove_multiples, generator=generator)


def elementary_vectors_from_matrix(M, remove_multiples=True, generator=False):
    r"""
    Compute elementary vectors of a subspace determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``remove_multiples`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:

    Return a list of elementary vectors in the kernel of ``M``.

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
        M = M.matrix_from_rows(M.pivot_rows()) # does not work for polynomial matrices
    except (ArithmeticError, NotImplementedError):
        warnings.warn('Could not determine rank of matrix. Expect wrong result!')

    rank, length = M.dimensions()
    minors = {}
    ring = M.base_ring()

    def ev_from_support(indices):
        r"""Return the elementary vector with support in a given list of indices."""
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
        evs = reduce_vectors_support(evs, generator=True)

    if generator:
        return evs
    return list(evs)


def elementary_vectors_from_minors(minors, dimensions, ring=None, remove_multiples=True, generator=False):
    r"""
    Compute elementary vectors determined by given maximal minors of a matrix.

    INPUT:

    - ``minors`` -- a list of maximal minors of a matrix

    - ``dim`` -- a tuple of the dimensions of the matrix corresponding to ``minors``

    - ``ring`` -- the ring in which the elementary vectors should live
                  by default, this is determined from ``minors``

    - ``remove_multiples`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix with
    maximal minors ``minors`` and dimension ``dimensions``.
    If ``generator`` is true, a generator of elementary vectors will be returned instead of a list.

    - If ``remove_multiples`` is true, the output will be reduced such that
      all elementary vectors have distinct support.
    
    - If ``generator`` is true, the output will be a generator object instead of a list.

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
    if not ring:
        ring = minors[0].base_ring()

    minors_dict = {}
    for k, indices in enumerate(Combinations(length, rank)):
        minors_dict[tuple(indices)] = minors[k]

    def ev_from_support(indices):
        r"""Return the elementary vector with support in a given list of indices."""
        element = zero_vector(ring, length)
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            minor = minors_dict[indices_minor]
            element[k] = (-1)**pos * minor
        return element

    evs = (ev_from_support(indices) for indices in Combinations(length, rank + 1))
    evs = (v for v in evs if v)
    if remove_multiples:
        evs = reduce_vectors_support(evs, generator=True)

    return evs if generator else list(evs)

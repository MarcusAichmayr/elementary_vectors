r"""Computing elementary vectors"""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import warnings
from sage.combinat.combination import Combinations
from sage.combinat.permutation import Permutation
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import zero_vector

from .utility import elementary_vector_from_indices, elementary_vector_from_indices_prevent_multiples, is_symbolic


def elementary_vectors(
    data,
    dimensions: tuple[int, int] = None,
    kernel: bool = True,
    generator: bool = False,
):
    r"""
    Compute elementary vectors of a subspace determined by a matrix or a list of maximal minors.

    INPUT:

    - ``data`` -- a matrix or a list of maximal minors

    - ``dimensions`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.

    - ``kernel`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix given by ``data``.
    The element ``data`` can be a matrix or a list of maximal minors.
    In the latter case, dimensions are also required.
    To compute the elementary vectors in the row space, pass false for ``kernel``.

    - If ``remove_multiples`` is true, the output will be reduced such that
      all elementary vectors have distinct support.

    - If ``generator`` is true, the output will be a generator object instead of a list.

    .. SEEALSO::

        :func:`~elementary_vectors_from_matrix`
        :func:`~elementary_vectors_from_minors`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors(M)
        [(0, 4, -2, -2, 0), (2, 0, 0, 0, -2)]
    
    By default, the elementary vectors in the kernel are computed.
    To consider the row space, pass ``kernel=False``::

        sage: elementary_vectors(M, kernel=False)
        [(0, -2, -4, 0, 0), (0, -2, 0, -4, 0), (-4, 0, 0, 0, -4), (0, 0, -2, 2, 0)]

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
            raise TypeError(
                "When computing elementary vectors from a list of "
                + "maximal minors, the dimensions of the corresponding "
                + "matrix are needed."
            )
        return elementary_vectors_from_minors(
            data,
            dimensions=dimensions,
            kernel=kernel,
            generator=generator,
        )
    return elementary_vectors_from_matrix(data, kernel=kernel, generator=generator)


def elementary_vectors_from_matrix(M, kernel: bool = True, generator: bool = False):
    r"""
    Compute elementary vectors of a subspace determined by the matrix ``M``.

    INPUT:

    - ``M`` -- a matrix

    - ``remove_multiples`` -- a boolean (default: ``True``)

    - ``kernel`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:

    Return a list of elementary vectors in the kernel of ``M``.
    To compute the elementary vectors in the row space, pass false for ``kernel``.

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
        sage: elementary_vectors_from_matrix(M, kernel=False)
        [(0, -2, -4, 0, 0), (0, -2, 0, -4, 0), (-4, 0, 0, 0, -4), (0, 0, -2, 2, 0)]
    """
    try:
        M = M.matrix_from_rows(M.pivot_rows())  # does not work for polynomial matrices
    except (ArithmeticError, NotImplementedError):
        warnings.warn("Could not determine rank of matrix. Expect wrong result!")

    rank, length = M.dimensions()
    minors = {}
    evs = (
        elementary_vector_from_indices_prevent_multiples(indices, minors, M, kernel=kernel)
        for indices in (Combinations(length, rank + 1) if kernel else Combinations(length, length - rank + 1))
    )
    evs = (v for v in evs if v)

    if generator:
        return evs
    return list(evs)


def elementary_vectors_from_minors(
    minors: list,
    dimensions: tuple[int, int],
    ring=None,
    kernel: bool = True,
    generator: bool = False,
):
    r"""
    Compute elementary vectors determined by given maximal minors of a matrix.

    INPUT:

    - ``minors`` -- a list of maximal minors of a matrix

    - ``dim`` -- a tuple of the dimensions of the matrix corresponding to ``minors``

    - ``ring`` -- the ring in which the elementary vectors should live
                  by default, this is determined from ``minors``

    - ``kernel`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix with
    maximal minors ``minors`` and dimension ``dimensions``.
    To compute the elementary vectors in the row space, pass false for ``kernel``.
    If ``generator`` is true, a generator of elementary vectors will be returned instead of a list.

    - If ``remove_multiples`` is true, the output will be reduced such that
      all elementary vectors have distinct support.

    - If ``generator`` is true, the output will be a generator object instead of a list.

    .. SEEALSO::

        :func:`~elementary_vectors`

    EXAMPLES::

        sage: from elementary_vectors.functions import elementary_vectors_from_minors
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]])
        sage: M
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
        sage: elementary_vectors_from_minors(m, [3, 5], kernel=False)
        [(0, -2, -4, 0, 0), (0, -2, 0, -4, 0), (-4, 0, 0, 0, -4), (0, 0, -2, 2, 0)]
    """
    rank, length = dimensions
    if not ring:
        ring = minors[0].base_ring()

    minors_dict = {}
    for k, indices in enumerate(Combinations(length, rank)):
        minors_dict[tuple(indices)] = minors[k]

    evs = (
        elementary_vector_from_indices_prevent_multiples(indices, minors_dict, length=length, ring=ring, kernel=kernel)
        for indices in (Combinations(length, rank + 1) if kernel else Combinations(length, length - rank + 1))
    )
    evs = (v for v in evs if v)

    return evs if generator else list(evs)


def kernel_matrix_using_elementary_vectors(M):
    """
    Division-free computation of the right kernel based on elementary vectors.

    INPUT:

    - ``M`` -- a matrix

    OUTPUT:
    A right kernel matrix of ``M``.
    Also works for symbolic matrices.

    .. NOTE::

        A ``ValueError`` is returned if ``M`` has no constant nonzero maximal minor.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, 1, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  1  1 -1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [1 0 0 1 1]
        [0 1 0 0 1]
        [0 0 1 1 2]
        sage: M = matrix([[1, 1, -1, -1, 0], [2, 1, -1, 0, 1], [1, 1, 1, 1, 1]])
        sage: M
        [ 1  1 -1 -1  0]
        [ 2  1 -1  0  1]
        [ 1  1  1  1  1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [-1  0  0 -1  2]
        [ 0 -1  1 -2  2]
        sage: var('a')
        a
        sage: M = matrix([[1, 0, 1, -1, 0], [0, 1, a, 1, -1]])
        sage: M
        [ 1  0  1 -1  0]
        [ 0  1  a  1 -1]
        sage: kernel_matrix_using_elementary_vectors(M)
        [    1     0     0     1     1]
        [    0     1     0     0     1]
        [    0     0     1     1 a + 1]

    TESTS::

        sage: kernel_matrix_using_elementary_vectors(identity_matrix(3, 3))
        []
        sage: _.dimensions()
        (0, 3)
        sage: kernel_matrix_using_elementary_vectors(matrix(0, 3))
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    try:
        M = M.matrix_from_rows(M.pivot_rows())  # does not work for polynomial matrices
    except (ArithmeticError, NotImplementedError):
        warnings.warn("Could not determine rank of matrix. Expect wrong result!")

    rank, length = M.dimensions()
    minors = {}

    if rank == length:
        return matrix(M.base_ring(), 0, length)

    constant_minor_found = False
    for indices_minor in Combinations(reversed(range(length)), rank):
        indices_minor = tuple(indices_minor)
        minors[indices_minor] = M.matrix_from_columns(indices_minor).det()
        minor = minors[indices_minor]
        if minor and not is_symbolic(minor):
            constant_minor_found = True
            break
    if not constant_minor_found:
        raise ValueError("Could not find a constant nonzero maximal minor.")

    return matrix(
        elementary_vector_from_indices(set(indices_minor).union([k]), minors, M)
        for k in range(length)
        if k not in indices_minor
    )


class EVs:
    r"""
    A class used to compute elementary vectors.
    
    Whenever a maximal minor is computed, it is stored in a dictionary for efficient reuse.
    Supports elementary vectors in the kernel and row space.

    EXAMPLES::

        sage: from elementary_vectors.functions import EVs
        sage: M = matrix([[1, 2, 3, 4, 5], [0, 1, 2, 2, 3]])
        sage: m = EVs(M)
        sage: m.elementary_vectors()
        [(1, -2, 1, 0, 0),
         (0, -2, 0, 1, 0),
         (1, -3, 0, 0, 1),
         (-2, 0, -2, 2, 0),
         (-1, 0, -3, 0, 2),
         (2, 0, 0, -3, 2),
         (0, -1, -1, 0, 1),
         (0, 0, 2, 1, -2)]
        sage: m.elementary_vectors(prevent_multiples=False)
        [(1, -2, 1, 0, 0),
         (0, -2, 0, 1, 0),
         (1, -3, 0, 0, 1),
         (-2, 0, -2, 2, 0),
         (-1, 0, -3, 0, 2),
         (2, 0, 0, -3, 2),
         (0, -2, 0, 1, 0),
         (0, -1, -1, 0, 1),
         (0, 2, 0, -1, 0),
         (0, 0, 2, 1, -2)]
        sage: m.elementary_vectors(kernel=False)
        [(-3, -1, 1, -2, 0), (2, 0, -2, 0, -2), (-2, -1, 0, -2, -1), (0, 1, 2, 2, 3)]
        sage: m.elementary_vectors(kernel=False, prevent_multiples=False)
        [(-3, -1, 1, -2, 0),
         (2, 0, -2, 0, -2),
         (-2, -1, 0, -2, -1),
         (1, 0, -1, 0, -1),
         (0, 1, 2, 2, 3)]
    """
    def __init__(self, M) -> None:
        try:
            self.matrix = M.matrix_from_rows(M.pivot_rows())  # does not work for polynomial matrices
        except (ArithmeticError, NotImplementedError):
            self.matrix = M
            warnings.warn("Could not determine rank of matrix. Expect wrong result!")

        self.length = self.matrix.ncols()
        self.rank = self.matrix.nrows()
        self.minors = {}
        self.minors_kernel = {}
        self.ring = self.matrix.base_ring()

        self.marked_minors = set()
        self.marked_minors_kernel = set()

    def minor_from_indices(self, indices, kernel: bool = True):
        r"""
        Compute a minor given by indices

        The minor is cached for efficient reuse.

        TESTS::

            sage: from elementary_vectors.functions import EVs
            sage: M = matrix([[1, 2, 3, 4, 5], [0, 1, 2, 2, 3]])
            sage: m = EVs(M)
            sage: m.minors
            {}
            sage: m.minors_kernel
            {}
            sage: m.minor_from_indices([0, 1])
            1
            sage: m.minors
            {(0, 1): 1}
            sage: m.minors_kernel
            {}
            sage: m.minor_from_indices([2, 3, 4], kernel=False)
            1
            sage: m.minors
            {(0, 1): 1}
            sage: m.minors_kernel
            {(2, 3, 4): 1}
            sage: m.minor_from_indices([0, 1, 3], kernel=False)
            1
            sage: m.minors
            {(0, 1): 1, (2, 4): -1}
            sage: m.minors_kernel
            {(0, 1, 3): 1, (2, 3, 4): 1}
            sage: m.minor_from_indices([2, 4])
            -1
            sage: m.minors
            {(0, 1): 1, (2, 4): -1}
            sage: m.minors_kernel
            {(0, 1, 3): 1, (2, 3, 4): 1}
        """
        indices = tuple(indices)
        try:
            if kernel:
                return self.minors[indices]
            return self.minor_of_kernel_matrix(indices)
        except KeyError:
            pass

        self.minors[indices] = self.matrix.matrix_from_columns(indices).det()
        return self.minors[indices]

    def minor_of_kernel_matrix(self, indices):
        r"""
        Return a minor of the kernel matrix given by indices

        Caches the result for efficient reuse and uses the corresponding minor
        """
        indices = tuple(indices)
        try:
            return self.minors_kernel[indices]
        except KeyError:
            pass
        indices_complement = tuple(i for i in range(self.length) if i not in indices)
        minor = self.minor_from_indices(indices_complement)
        self.minors_kernel[indices] = Permutation(i + 1 for i in list(indices_complement) + list(indices)).sign() * minor
        return self.minors_kernel[indices]

    def elementary_vector(self, indices: list, kernel: bool = True):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        OUTPUT:
        If ``allow_multiple`` is false, a ``ValueError`` will be raised if a zero minor is reused.
        """
        element = zero_vector(self.ring, self.length)
        nonzero_detected = False
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            minor = self.minor_from_indices(indices_minor, kernel=kernel)
            if minor == 0:
                continue
            nonzero_detected = True
            element[k] = (-1) ** pos * minor
        if nonzero_detected:
            return element
        raise ValueError("Indices correspond to zero vector!")

    def elementary_vector_prevent_multiple(self, indices: list, kernel: bool = True):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        OUTPUT:
        If ``allow_multiple`` is false, a ``ValueError`` will be raised if a zero minor is reused.
        """
        element = zero_vector(self.ring, self.length)
        nonzero_detected = False
        zero_minors = []
        multiple_detected = False
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            if kernel:
                if indices_minor in self.marked_minors_kernel:
                    multiple_detected = True
            else:
                if indices_minor in self.marked_minors:
                    multiple_detected = True
            minor = self.minor_from_indices(indices_minor, kernel=kernel)
            if minor == 0:
                zero_minors.append(indices_minor)
                continue
            nonzero_detected = True
            element[k] = (-1) ** pos * minor
        if nonzero_detected:
            for marked_minor in zero_minors:
                if kernel:
                    self.marked_minors_kernel.add(marked_minor)
                else:
                    self.marked_minors.add(marked_minor)
                if multiple_detected:
                    raise ValueError("Multiple detected!")
            return element
        raise ValueError("Indices correspond to zero vector!")

    def reset_set_for_preventing_multiples(self, kernel: bool = True) -> None:
        if kernel:
            self.marked_minors_kernel = set()
        else:
            self.marked_minors = set()

    def elementary_vectors_generator(self, kernel: bool = True, prevent_multiples: bool = True) -> list:
        if prevent_multiples:
            self.reset_set_for_preventing_multiples(kernel=kernel)
        for indices in (Combinations(self.length, self.rank + 1) if kernel else Combinations(self.length, self.length - self.rank + 1)):
            try:
                if prevent_multiples:
                    yield self.elementary_vector_prevent_multiple(indices, kernel=kernel)
                else:
                    yield self.elementary_vector(indices, kernel=kernel)
            except ValueError:
                pass

    def elementary_vectors(self, kernel: bool = True, prevent_multiples: bool = True) -> list:
        return list(self.elementary_vectors_generator(kernel=kernel, prevent_multiples=prevent_multiples))
        # result = []
        # if prevent_multiples:
        #     self.reset_set_for_preventing_multiples(kernel=kernel)
        # for indices in (Combinations(self.length, self.rank + 1) if kernel else Combinations(self.length, self.length - self.rank + 1)):
        #     try:
        #         if prevent_multiples:
        #             result.append(self.elementary_vector_prevent_multiple(indices, kernel=kernel))
        #         else:
        #             result.append(self.elementary_vector(indices, kernel=kernel))
        #     except ValueError:
        #         pass
        # return result

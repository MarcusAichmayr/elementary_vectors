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
from collections.abc import Generator
from sage.combinat.combination import Combinations
from sage.combinat.permutation import Permutation
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import zero_vector

from .utility import is_symbolic


def elementary_vectors(M, kernel: bool = True, generator: bool = False):
    r"""
    Compute elementary vectors of a subspace determined by a matrix or a list of maximal minors.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    - ``generator`` -- a boolean (default: ``False``)

    OUTPUT:
    Return a list of elementary vectors in the kernel of a matrix given by the matrix ``M``.
    To compute the elementary vectors in the row space, pass false for ``kernel``.

    - If ``generator`` is true, the output will be a generator object instead of a list.

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]])
        sage: M
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

    TESTS::

        sage: elementary_vectors(random_matrix(QQ, 0, 4))
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        sage: elementary_vectors(matrix([[1, 0, 0], [0, 0, 0]]))
        [(0, -1, 0), (0, 0, -1)]
    """
    if generator:
        return ElementaryVectors(M).generator(kernel=kernel)
    return ElementaryVectors(M).elements(kernel=kernel)


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
    rank, length = M.dimensions()

    if rank == length:
        return matrix(M.base_ring(), 0, length)

    evs = ElementaryVectors(M)
    constant_minor_found = False
    for indices in Combinations(reversed(range(length)), rank):
        minor = evs.minor(indices)
        if minor and not is_symbolic(minor):
            constant_minor_found = True
            break
    if not constant_minor_found:
        raise ValueError("Could not find a constant nonzero maximal minor.")

    return matrix(evs.element(set(indices).union([k])) for k in range(length) if k not in indices)


class ElementaryVectors:
    r"""
    A class used to compute elementary vectors.
    
    Whenever a maximal minor is computed, it is stored in a dictionary for efficient reuse.
    Supports elementary vectors in the kernel and row space.

    EXAMPLES::

        sage: from elementary_vectors.functions import ElementaryVectors
        sage: M = matrix([[1, 2, 3, 4, 5], [0, 1, 2, 2, 3]])
        sage: evs = ElementaryVectors(M)
        sage: evs.elements()
        [(1, -2, 1, 0, 0),
         (0, -2, 0, 1, 0),
         (1, -3, 0, 0, 1),
         (-2, 0, -2, 2, 0),
         (-1, 0, -3, 0, 2),
         (2, 0, 0, -3, 2),
         (0, -1, -1, 0, 1),
         (0, 0, 2, 1, -2)]
        sage: evs.elements(prevent_multiples=False)
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
        sage: evs.elements(kernel=False)
        [(-3, -1, 1, -2, 0), (2, 0, -2, 0, -2), (-2, -1, 0, -2, -1), (0, 1, 2, 2, 3)]
        sage: evs.elements(kernel=False, prevent_multiples=False)
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
        self.ring = self.matrix.base_ring()
        self.minors = {}
        self.minors_dual = {}
        self.set_combinations()
        self.set_combinations_dual()

        self.marked_minors = set()
        self.marked_minors_dual = set()

    def minor(self, indices, kernel: bool = True):
        r"""
        Compute a minor given by indices

        The minor is cached for efficient reuse.

        TESTS::

            sage: from elementary_vectors.functions import ElementaryVectors
            sage: M = matrix([[1, 2, 3, 4, 5], [0, 1, 2, 2, 3]])
            sage: evs = ElementaryVectors(M)
            sage: evs.minors
            {}
            sage: evs.minors_dual
            {}
            sage: evs.minor([0, 1])
            1
            sage: evs.minors
            {(0, 1): 1}
            sage: evs.minors_dual
            {}
            sage: evs.minor([2, 3, 4], kernel=False)
            1
            sage: evs.minors
            {(0, 1): 1}
            sage: evs.minors_dual
            {(2, 3, 4): 1}
            sage: evs.minor([0, 1, 3], kernel=False)
            1
            sage: evs.minors
            {(0, 1): 1, (2, 4): -1}
            sage: evs.minors_dual
            {(0, 1, 3): 1, (2, 3, 4): 1}
            sage: evs.minor([2, 4])
            -1
            sage: evs.minors
            {(0, 1): 1, (2, 4): -1}
            sage: evs.minors_dual
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
            return self.minors_dual[indices]
        except KeyError:
            pass
        indices_dual = tuple(i for i in range(self.length) if i not in indices)
        minor = self.minor(indices_dual)
        self.minors_dual[indices] = Permutation(i + 1 for i in list(indices_dual) + list(indices)).sign() * minor
        return self.minors_dual[indices]

    def set_combinations(self, combinations=None) -> None:
        r"""Set or reset combinations."""
        if combinations is None:
            self._combinations = Combinations(self.length, self.rank + 1)
        else:
            self._combinations = combinations

    def set_combinations_dual(self, combinations=None) -> None:
        r"""Set or reset combinations."""
        if combinations is None:
            self._combinations_dual = Combinations(self.length, self.length - self.rank + 1)
        else:
            self._combinations_dual = combinations

    def element(self, indices: list, kernel: bool = True):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            Raises a ``ValueError`` if the indices correspond to the zero vector.
        """
        element = self._zero_element()
        nonzero_detected = False
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            minor = self.minor(indices_minor, kernel=kernel)
            if minor == 0:
                continue
            nonzero_detected = True
            element[k] = (-1) ** pos * minor
        if nonzero_detected:
            return element
        raise ValueError("Indices correspond to zero vector!")

    def element_prevent_multiple(self, indices: list, kernel: bool = True):
        r"""
        Compute the elementary vector corresponding to a list of indices.

        .. NOTE::

            If this results in a multiple of a previous element, a ``ValueError`` is raised.
        """
        element = self._zero_element()
        nonzero_detected = False
        zero_minors = []
        multiple_detected = False
        for pos, k in enumerate(indices):
            indices_minor = tuple(i for i in indices if i != k)
            if kernel:
                if indices_minor in self.marked_minors:
                    multiple_detected = True
            else:
                if indices_minor in self.marked_minors_dual:
                    multiple_detected = True
            minor = self.minor(indices_minor, kernel=kernel)
            if minor == 0:
                zero_minors.append(indices_minor)
                continue
            nonzero_detected = True
            element[k] = (-1) ** pos * minor
        if nonzero_detected:
            for marked_minor in zero_minors:
                if kernel:
                    self.marked_minors.add(marked_minor)
                else:
                    self.marked_minors_dual.add(marked_minor)
                if multiple_detected:
                    raise ValueError("Multiple detected!")
            return element
        raise ValueError("Indices correspond to zero vector!")

    def random_element(self, kernel: bool = True):
        r"""
        Return a random elementary vector

        .. NOTE::

            If no elementary vector exists or the zero vector has been generated, ``None`` is returned.
        """
        try:
            if kernel:
                return self.element(self._combinations.random_element())
            return self.element(self._combinations_dual.random_element(), kernel=False)
        except ValueError: # no elementary vectors exist or generated zero vector
            return

    def _zero_element(self):
        return zero_vector(self.ring, self.length)

    def _reset_set_for_preventing_multiples(self, kernel: bool = True) -> None:
        if kernel:
            self.marked_minors = set()
        else:
            self.marked_minors_dual = set()

    def generator(
        self,
        kernel: bool = True,
        prevent_multiples: bool = True,
        reverse: bool = False
    ) -> Generator:
        r"""Return a generator of elementary vectors"""
        if prevent_multiples:
            self._reset_set_for_preventing_multiples(kernel=kernel)
        combinations = self._combinations if kernel else self._combinations_dual
        if reverse:
            combinations = reversed(combinations)
        for indices in combinations:
            try:
                if prevent_multiples:
                    yield self.element_prevent_multiple(indices, kernel=kernel)
                else:
                    yield self.element(indices, kernel=kernel)
            except ValueError:
                pass

    def elements(self, kernel: bool = True, prevent_multiples: bool = True) -> list:
        r"""Return a list of elementary vectors"""
        return list(self.generator(kernel=kernel, prevent_multiples=prevent_multiples))

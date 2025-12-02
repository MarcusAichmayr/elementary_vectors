r"""
Supports of circuits and cocircuits
===================================

For some applications, it is sufficient to compute only the supports of the circuits.
Compare also with the circuits of linear matroids.
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

from typing import List

from sage.matrix.constructor import Matrix

from .elements import CircuitEnumerator


def circuit_supports(matrix: Matrix) -> List[List[int]]:
    r"""
    Compute the supports of the circuits of a matrix.

    .. SEEALSO::

        :func:`cocircuit_supports`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: circuit_supports(M)
        [[0, 1, 2], [0, 1, 3], [2, 3]]

    TESTS:

    This generates the empty list::

        sage: M = matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
        sage: M
        [1 0 0 0]
        [0 1 0 0]
        sage: circuit_supports(M)
        [[2], [3]]
    """
    return CircuitSupportEnumerator(matrix).circuits()


def cocircuit_supports(matrix: Matrix) -> List[List[int]]:
    r"""
    Compute the supports of the cocircuits of a matrix.

    .. SEEALSO::

        :func:`circuit_supports`

    EXAMPLES::

        sage: from elementary_vectors import *
        sage: M = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
        sage: M
        [1 2 0 0]
        [0 1 2 3]
        sage: cocircuit_supports(M)
        [[1, 2, 3], [0, 2, 3], [0, 1]]
    """
    return CircuitSupportEnumerator(matrix).cocircuits()


class CircuitSupportEnumerator(CircuitEnumerator):
    r"""
    Circuit support enumerator.
    """
    def _compute_minor(self, indices: List[int]) -> int:
        minor = super()._compute_minor(indices)
        return 0 if minor == 0 else 1

    def _repr_(self) -> str:
        return f"Circuit support enumerator of {self.rank}x{self.length} matrix"

    def _zero_element(self) -> List[int]:
        return []

    @staticmethod
    def _set_element_entry(element: List[int], index: int, value: int) -> None:
        element.append(index)

    @staticmethod
    def _is_element_zero(element: List[int]) -> bool:
        return len(element) == 0

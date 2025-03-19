"""Utility functions"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from collections.abc import Generator

from sage.combinat.combination import Combinations
from sage.functions.generalized import sign
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.structure.sage_object import SageObject

from elementary_vectors.functions import ElementaryVectors


def solve_left_for_roots(A, b):
    r"""
    Find a solution for ``x*A = b`` that works for matrices with roots.

    INPUT:

    - ``A`` -- a matrix

    - ``b`` -- a vector

    NOTE::

        The built in method ``solve_left`` for matrices fails occasionally.
    """
    M = matrix(list(A) + [-b]).T.right_kernel_matrix()
    x = matrix(M.column(-1)).solve_right(vector([1]))
    return (x * M)[:-1]


def solve_without_division(A, b):
    r"""
    Solve a linear system of equations without division.

    The system is ``A x = c b`` where ``c`` is a positive constant.
    Uses an elementary vector.

    EXAMPLES::

        sage: from vectors_in_intervals.utility import solve_without_division
        sage: A = matrix([[1, 2], [0, 1], [1, -1]])
        sage: b = vector([1, 0, 1])
        sage: solve_without_division(A, b)
        (1, 0)
        sage: A = matrix([[1, 4], [0, 2], [1, -2]])
        sage: b = vector([6, 2, 0])
        sage: solve_without_division(A, b)
        (4, 2)
        sage: A.solve_right(b)
        (2, 1)
        sage: A = matrix([[1, 1, 1], [0, 1, 2]])
        sage: b = vector([2, 3])
        sage: solve_without_division(A, b)
        (0, 1, 1)
    """
    Ab = matrix.block([[A, matrix.column(b)]])
    ev = next(ElementaryVectors(Ab).generator(reverse=True))
    return -sign(ev[-1]) * ev[:-1]


class CombinationsIncluding(SageObject):
    r"""
    Combinatorial object of all combinations that include given elements

    EXAMPLES:

    We generate all subsets of ``range(4)`` with ``2`` elements that include the element ``2``::

        sage: from vectors_in_intervals.utility import CombinationsIncluding
        sage: C = CombinationsIncluding(4, 2, [2])
        sage: list(C)
        [[0, 2], [1, 2], [2, 3]]
        sage: list(reversed(C))
        [[2, 3], [1, 2], [0, 2]]
    """
    def __init__(self, mset, k, elements=None):
        self.combinations = Combinations(set(range(mset)) - set(elements), k - len(elements))
        if elements is None:
            self.elements = []
        else:
            self.elements = elements

    def __iter__(self) -> Generator[list]:
        for combination in self.combinations:
            yield sorted(list(self.elements) + combination)
    
    def __reversed__(self) -> Generator[list]:
        for combination in reversed(self.combinations):
            yield sorted(list(self.elements) + combination)

    def random_element(self) -> list:
        return sorted(list(self.elements) + self.combinations.random_element())

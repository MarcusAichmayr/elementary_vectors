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

from sage.modules.free_module_element import vector


def kernel_vector_support_given(M, indices: list):
    r"""
    Return a right kernel vector such that the support is a subset of given indices.

    INPUT:

    - ``M`` -- a matrix

    - ``indices`` -- a list of indices

    OUTPUT:
    a vector in the right kernel of ``M`` such that the support is a subset of ``indices``.

    EXAMPLES::

        sage: from elementary_vectors.utility import kernel_vector_support_given
        sage: M = matrix([[1, 2, 0, 0], [0, 1, -1, 0]])
        sage: v = kernel_vector_support_given(M, [0, 1, 2])
        sage: max(v, -v)
        (2, -1, -1, 0)
        sage: kernel_vector_support_given(M, [3])
        (0, 0, 0, 1)
        sage: kernel_vector_support_given(M, [0, 3])
        (0, 0, 0, 1)
    """
    try:
        subvector_list = list(M.matrix_from_columns(indices).right_kernel_matrix()[0])
    except IndexError as exc:
        raise ValueError(
            "Right kernel restricted to column ``indices`` is empty."
        ) from exc
    for k in range(M.ncols()):
        if not k in indices:
            subvector_list.insert(k, 0)
    return vector(subvector_list)


def is_symbolic(expression):
    r"""
    Return whether this element is a symbolic expression.

    If it belongs to the symbolic ring but doesn't contain any variables it does not count as "symbolic".

    EXAMPLES::

        sage: from elementary_vectors.utility import is_symbolic
        sage: is_symbolic(5)
        False
        sage: var('a, b')
        (a, b)
        sage: is_symbolic(a)
        True
        sage: is_symbolic(-a)
        True
        sage: is_symbolic(b^2 - a)
        True
        sage: is_symbolic(SR(5))
        False
    """
    if hasattr(expression, "variables"):
        return bool(expression.variables())
    return False

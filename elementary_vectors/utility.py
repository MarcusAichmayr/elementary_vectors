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

from sign_vectors import sign_vector


def kernel_vector_support_given(M, indices):
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


def is_symbolic(value):
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
    if hasattr(value, "variables"):
        return bool(value.variables())
    return False


def conformal_elimination(vector1, vector2, indices: list[int] = None):
    r"""
    Apply conformal elimination to two real vectors to find a new vector.

    INPUT:

    - ``vector1`` -- a real vector

    - ``vector2`` -- a real vector

    - ``indices`` -- a list of indices (default: ``None``)

    OUTPUT:
    Return a new vector ``z = x + a y`` where ``a > 0``, such that ``z[e] == 0``
    for some ``e`` in ``indices`` and ``Z_k <= X_k`` for ``k`` in ``indices``
    and ``Z_f = (X o Y)_f`` for ``f`` not in ``D(X, Y)``.
    Here, ``X``, ``Y`` and ``Z`` are the sign vectors
    corresponding to ``x``, ``y`` and ``z``.

    .. NOTE::

        If ``indices`` is not given, the whole list of separating elements
        will be considered instead. (default)

    EXAMPLES::

        sage: from elementary_vectors.utility import conformal_elimination
        sage: x = vector([1, 0, 2])
        sage: y = vector([-1, 1, 1])
        sage: conformal_elimination(x, y)
        (0, 1, 3)
    """
    if indices is None:
        indices = []
    if vector1.length() != vector2.length():
        raise ValueError("Vectors have different length.")
    separating_elements = sign_vector(vector1).separating_elements(sign_vector(vector2))
    if not separating_elements:
        raise ValueError("List of separating elements is empty.")
    if not indices:
        indices = separating_elements
    elif not all(s in separating_elements for s in indices):
        raise ValueError("indices is not a subset of separating_elements.")
    lam = max(
        vector1[e] / vector2[e] for e in indices
    )  # lam < 0 since e in separating_elements
    return vector1 - lam * vector2

r"""Functions for working with oriented matroids"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.posets.posets import Poset

from sign_vectors import SignVector, sign_vector, zero_sign_vector

from .utility import exclude_indices


def closure(iterable) -> set[SignVector]:
    r"""
    Compute the closure of given sign vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

    OUTPUT:
    Return the closure of ``iterable``.

    .. NOTE::

       The sign vector :math:`X` is in the closure
       of a set of sign vectors :math:`W`
       if there exists :math:`Y \in W` with :math:`X \leq Y`.

    EXAMPLES:

    We consider a list consisting of only one sign vector::

        sage: from sign_vectors import *
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: closure(W)
        {(000), (+00), (0-0), (+-0)}

    Now, we consider a list of three sign vectors::

        sage: W = [sign_vector("++-"), sign_vector("-00"), sign_vector("0--")]
        sage: W
        [(++-), (-00), (0--)]
        sage: closure(W)
        {(000), (-00), (00-), (0+0), (0+-), (+00), (+0-), (++0), (++-), (0-0), (0--)}

    TESTS::

        sage: closure([])
        set()
    """
    if not iterable:
        return set()
    for _ in iterable:
        length = _.length()
        break

    output = [{zero_sign_vector(length)}]
    new_elements = set()
    for i in range(length):
        X = sign_vector(1 if k == i else 0 for k in range(length))
        for Z in iterable:
            if X <= Z:
                new_elements.add(X)
                break
        Y = sign_vector(-1 if k == i else 0 for k in range(length))
        for Z in iterable:
            if Y <= Z:
                new_elements.add(Y)
                break
    output.append(new_elements)
    for i in range(1, length + 1):
        new_elements = set()
        for X in output[1]:  # X has always |supp(X)| = 1
            for Y in output[i]:
                # TODO: utilize that the supports are sorted
                if len(set(X.support() + Y.support())) == i + 1:
                    Z = X.compose(Y)
                    if Z not in new_elements and any(Z <= V for V in iterable):
                        new_elements.add(Z)
        if not new_elements:
            break
        output.append(new_elements)
        if len(new_elements) == 1:
            break
    return set().union(*output)


def contraction(iterable, indices: list[int], keep_components: bool = False) -> set:
    r"""
    Return all sign vectors or vectors that are zero on given components.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors

    - ``indices`` -- a list of indices.

    - ``keep_components`` -- a boolean

    OUTPUT:

    - If ``keep_components`` is false, remove entries in ``indices``. (default)

    - If ``keep_components`` is true, keep entries in ``indices``.

    EXAMPLES::

        sage: from sign_vectors import *
        sage: W = [sign_vector("++0"), sign_vector("-00"), sign_vector("00+")]
        sage: W
        [(++0), (-00), (00+)]

    Only the third sign vector has a zero at the component with index ``0``.
    Removing this component leads to the following result::

        sage: contraction(W, [0])
        {(0+)}
        sage: contraction(W, [1])
        {(-0), (0+)}
        sage: contraction(W, [2])
        {(-0), (++)}

    The second sign vector has zeros at positions ``1`` and ``2``::

        sage: contraction(W, [1, 2])
        {(-)}

    We take the examples from before. With ``keep_components=True``, we keep the
    zero components of the appropriate sign vectors::

        sage: contraction(W, [0], keep_components=True)
        {(00+)}
        sage: contraction(W, [1], keep_components=True)
        {(-00), (00+)}
        sage: contraction(W, [2], keep_components=True)
        {(-00), (++0)}
        sage: contraction(W, [1, 2], keep_components=True)
        {(-00)}

    This function also works for matrices or lists of vectors::

        sage: l = [vector([0, 0, 1]), vector([0, 2, 1]), vector([-1, 0, 1])]
        sage: contraction(l, [0])
        {(0, 1), (2, 1)}
        sage: A = matrix([[1, 1, 0], [0, 1, 0]])
        sage: contraction(A, [2])
        {(0, 1), (1, 1)}
    """
    if not iterable:
        return iterable
    if keep_components:

        def vec(iterable):
            return iterable
    else:
        vec = exclude_indices(iterable, indices)

    return set(vec(X) for X in iterable if not any(e in indices for e in X.support()))


def deletion(iterable, indices: list[int]) -> set:
    r"""
    Remove given components from an iterable of sign vectors or vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors

    - ``indices`` -- a list of indices

    EXAMPLES::

        sage: from sign_vectors import *
        sage: W = [sign_vector("++0"), sign_vector("00-"), sign_vector("+00")]
        sage: W
        [(++0), (00-), (+00)]
        sage: deletion(W, [0])
        {(00), (+0), (0-)}

    Duplicate sign vectors are removed if they would occur::

        sage: deletion(W, [1])
        {(+0), (0-)}
        sage: deletion(W, [1, 2])
        {(0), (+)}

    This function also works for lists of vectors::

        sage: l = [vector([0, 0, 1]), vector([0, 2, 1]), vector([-1, 0, 1])]
        sage: deletion(l, [1])
        {(-1, 1), (0, 1)}
    """
    if not iterable:
        return iterable

    return set(exclude_indices(iterable, indices)(X) for X in iterable)


def plot_sign_vectors(iterable, vertex_size: int = 600, figsize: int = None, aspect_ratio=None):
    r"""
    Plot the Hasse Diagram of sign vectors using the conformal relation.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

    - ``vertex_size`` -- the size of the vertices in the plot (default: 600)

    - ``figsize`` -- the size of the figure (default: None)

    - ``aspect_ratio`` -- the aspect ratio of the plot (default: None)
    """
    Poset((iterable, lambda X, Y: X.conforms(Y))).plot(
        vertex_size=vertex_size,
        element_color="white",
        vertex_shape="",
    ).show(figsize=figsize, aspect_ratio=aspect_ratio)

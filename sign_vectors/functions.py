r"""Functions for working with oriented matroids."""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.misc.flatten import flatten
from sage.combinat.posets.posets import Poset

from .utility import _subvector
from sign_vectors import sign_vector, zero_sign_vector


def closure(W, separate=False):
    r"""
    Compute the closure of a list of sign vectors.

    INPUT:

    - ``W`` -- an iterable of sign vectors

    - ``separate`` -- boolean (default: ``False``)

    OUTPUT:

    If ``separate`` is false, return the closure of ``W``. (default)

    If ``separate`` is true, separate the closure into sets, where each element
    has the same number of zero entries.

    .. NOTE::

       The sign vector :math:`X` is in the closure
       of a set of sign vectors :math:`W`
       if there exists :math:`Y \in W` with :math:`X \leq Y`.

    EXAMPLES:

    We consider a list consisting of only one sign vector::

        sage: from sign_vectors import sign_vector, closure
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: closure(W)
        {(000), (+00), (0-0), (+-0)}

    With the optional argument ``separate=True``, we can separate the resulting
    list into three sets.
    Each sign vector in such a set has the same number of zero entries::

        sage: closure(W, separate=True)
        [{(000)}, {(+00), (0-0)}, {(+-0)}]

    Now, we consider a list of three sign vectors::

        sage: W = [sign_vector("++-"), sign_vector("-00"), sign_vector("0--")]
        sage: W
        [(++-), (-00), (0--)]
        sage: closure(W)
        {(000), (-00), (00-), (+00), (+0-), (++0), (++-), (0--), (0+0), (0+-), (0-0)}
        sage: closure(W, separate=True)
        [{(000)},
         {(-00), (00-), (+00), (0+0), (0-0)},
         {(0--), (0+-), (+0-), (++0)},
         {(++-)}]

    TESTS::

        sage: closure([])
        set()
    """
    if not W:
        return set()
    for _ in W:
        n = _.length()
        break

    F = [{zero_sign_vector(n)}]
    F_new = set()
    for i in range(n):
        X = sign_vector(1 if k == i else 0 for k in range(n))
        for Z in W:
            if X <= Z:
                F_new.add(X)
                break
        Y = sign_vector(-1 if k == i else 0 for k in range(n))
        for Z in W:
            if Y <= Z:
                F_new.add(Y)
                break
    F.append(F_new)
    for i in range(1, n + 1):
        F_new = set()
        for X in F[1]:  # X has always |supp(X)| = 1
            for Y in F[i]:
                if len(set(X.support() + Y.support())) == i + 1:  # TODO: utilize that the supports are sorted
                    Z = X.compose(Y)
                    if Z not in F_new:
                        if any(Z <= V for V in W):
                            F_new.add(Z)
        if F_new == set():
            break
        else:
            F.append(F_new)
            if len(F_new) == 1:
                break
    if separate:
        return F
    else:
        return set().union(*F)


def contraction(F, R, keep_components=False):
    r"""
    Return all sign vectors that are zero on ``R``. Also works for real vectors.

    INPUT:

    - ``F`` -- an iterable of sign vectors or vectors; or a matrix.

    - ``R`` -- a list of indices.

    - ``keep_components`` -- a boolean (default: ``False``).

    OUTPUT:

    - If ``keep_components`` is false, remove entries in ``R``. (default)

    - If ``keep_components`` is true, keep entries in ``R``.

    EXAMPLES::

        sage: from sign_vectors import sign_vector, contraction
        sage: W = [sign_vector("++0"), sign_vector("-00"), sign_vector("00+")]
        sage: W
        [(++0), (-00), (00+)]

    Only the third sign vector has a zero at the component with index ``0``.
    Removing this component leads to the following result::

        sage: contraction(W, [0])
        {(0+)}
        sage: contraction(W, [1])
        {(0+), (-0)}
        sage: contraction(W, [2])
        {(++), (-0)}

    The second sign vector has zeros at positions ``1`` and ``2``::

        sage: contraction(W, [1, 2])
        {(-)}

    We take the examples from before. With ``keep_components=True``, we keep the
    zero components of the appropriate sign vectors::

        sage: contraction(W, [0], keep_components=True)
        {(00+)}
        sage: contraction(W, [1], keep_components=True)
        {(00+), (-00)}
        sage: contraction(W, [2], keep_components=True)
        {(++0), (-00)}
        sage: contraction(W, [1, 2], keep_components=True)
        {(-00)}

    This function also works for matrices or lists of vectors::

        sage: l = [vector([0,0,1]), vector([0,2,1]), vector([-1,0,1])]
        sage: contraction(l, [0])
        {(0, 1), (2, 1)}
        sage: A = matrix([[1,1,0],[0,1,0]])
        sage: contraction(A, [2])
        {(0, 1), (1, 1)}
    """
    if F == []:
        return F

    if keep_components:
        def vec(v):
            return v
    else:
        vec = _subvector(F, R)

    return set(vec(X) for X in F if not any(e in R for e in X.support()))


def deletion(F, R):
    r"""
    Remove the components corresponding to ``R`` from a list of sign vectors.
    Also works for real vectors.

    INPUT:

    - ``F`` -- an iterable of sign vectors or real vectors; or a matrix.

    - ``R`` -- a list of indices.

    EXAMPLES::

        sage: from sign_vectors import sign_vector, deletion
        sage: W = [sign_vector("++0"), sign_vector("00-"), sign_vector("+00")]
        sage: W
        [(++0), (00-), (+00)]
        sage: deletion(W, [0])
        {(00), (0-), (+0)}

    Duplicate sign vectors are removed if they would occur::

        sage: deletion(W, [1])
        {(0-), (+0)}
        sage: deletion(W, [1, 2])
        {(0), (+)}

    This function also works for lists of vectors::

        sage: l = [vector([0,0,1]), vector([0,2,1]), vector([-1,0,1])]
        sage: deletion(l, [1])
        {(-1, 1), (0, 1)}
    """
    if F == []:
        return F

    return set(_subvector(F, R)(X) for X in F)


def plot_sign_vectors(L, vertex_size=600, figsize=10, aspect_ratio=4/8):
    r"""Plot the Hasse Diagram of a list of given sign vectors using the conformal relation."""
    fcn = lambda X, Y: X.conforms(Y)
    P = Poset((L, fcn))
    P.plot(vertex_size=vertex_size, element_color='white', cover_style='-', vertex_shape='+').show(figsize=figsize, aspect_ratio=aspect_ratio)

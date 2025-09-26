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

from . import SignVector, zero_sign_vector


def lower_closure(iterable: set[SignVector]) -> set[SignVector]:
    r"""
    Compute the lower closure of given sign vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

    OUTPUT:
    Return the lower closure of ``iterable`` as a set of sign vectors.

    .. NOTE::

       The sign vector :math:`X` is in the lower closure
       of a set of sign vectors :math:`W`
       if there exists :math:`Y \in W` with :math:`X \leq Y`.

    EXAMPLES:

    We consider a list consisting of only one sign vector::

        sage: from sign_vectors import *
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: lower_closure(W)
        {(000), (+00), (0-0), (+-0)}

    Now, we consider a list of three sign vectors::

        sage: W = [sign_vector("++-"), sign_vector("-00"), sign_vector("0--")]
        sage: W
        [(++-), (-00), (0--)]
        sage: lower_closure(W)
        {(000), (-00), (00-), (0+0), (0+-), (+00), (+0-), (++0), (++-), (0-0), (0--)}

    TESTS::

        sage: lower_closure([])
        set()
    """
    if not iterable:
        return set()

    for _ in iterable:
        length = _.length()
        break

    max_support_length = 0
    same_support_list = [{zero_sign_vector(length)}]

    for sv in iterable:
        support_length = len(sv.support())

        while support_length > max_support_length:
            same_support_list.append(set())
            max_support_length += 1

        same_support_list[support_length].add(sv)

    for i in range(max_support_length, 1, -1):
        for sv in same_support_list[i]:
            for s in sv.support():
                same_support_list[i - 1].add(sv.set_to_zero([s]))

    return set().union(*same_support_list)

def upper_closure(iterable) -> tuple[set[SignVector], list[int]]:
    r"""
    Compute the upper closure of given sign vectors.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

    OUTPUT:
    Return the upper closure of ``iterable`` as a set of relevant sign vectors 
    and a list of irrelevant components. Note, that these can be unpacked with
    the function unpack_irrelevant_components.

    .. NOTE::

       The sign vector :math:`X` is in the upper closure
       of a set of sign vectors :math:`W`
       if there exists :math:`Y \in W` with :math:`X \geq Y`.

    EXAMPLES:

    We consider a list consisting of only one sign vector::

        sage: from sign_vectors import *
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: upper_closure(W)
        ({(+-0)}, [2])

    Now, we consider a list of three sign vectors::

        sage: W = [sign_vector("++-"), sign_vector("-00"), sign_vector("0--")]
        sage: W
        [(++-), (-00), (0--)]
        sage: upper_closure(W)
        ({(-00), (---), (++-), (-+0), (--+), (-++), (--0), (0--), (-0-), (+--), (-+-), (-0+)}, [])

    TESTS::

        sage: lower_closure([])
        set()
    """
    if not iterable:
        return set()

    for _ in iterable:
        length = _.length()
        irrelevant_c = _._zero_support()
        break

    min_support_length = length
    same_support_list = [set()]

    for sv in iterable:
        support_length = len(sv.support())
        irrelevant_c = irrelevant_c & sv._zero_support()

        while support_length < min_support_length:
            same_support_list.append(set())
            min_support_length -= 1

        same_support_list[length-support_length].add(sv)

    irrelevant_len = len(list(irrelevant_c))

    for i in range(length - min_support_length, irrelevant_len, -1):
        for sv in same_support_list[i]:
            for s in list(sv._zero_support()-irrelevant_c):
                same_support_list[i - 1].add(sv.set_to_plus([s]))
                same_support_list[i - 1].add(sv.set_to_minus([s]))

    return set().union(*same_support_list), list(irrelevant_c)

def unpack_irrelevant_components(iterable, irrelevant_components) -> list[SignVector]:
    r"""
    Compute the upper closure of given sign vectors.

    INPUT:

    - ``iterable`` -- a set of sign vectors
    - ``irrelevant_components`` -- list of indices

    OUTPUT:
    Return the full set of ``iterable``, with unpacked irrelevant components.

    EXAMPLES:

    We consider an output of the upper_closure function::

        sage: from sign_vectors import *
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: U = upper_closure(W); U
        ({(+-0)}, [2])
        sage: unpack_irrelevant_components(*U)
        [(+-0), (+-+), (+--)]

    Now, we consider a larger example of irrelevant components::

        sage: from sign_vectors import *
        sage: W = [sign_vector("+0-00"), sign_vector("+0000"), sign_vector("+0++0")]
        sage: W
        [(+0-00), (+0000), (+0++0)]
        sage: irrelevant_components = [1,4]
        sage: unpack_irrelevant_components(W, irrelevant_components)
        [(+0000), (+0++0), (+0-00), (++000), (+-000), (++++0), (+-++0), (++-00),
            (+--00), (+000+), (+000-), (+0+++), (+0++-), (+0-0+), (+0-0-), (++00+),
            (++00-), (+-00+), (+-00-), (+++++), (++++-), (+-+++), (+-++-), (++-0+),
            (++-0-), (+--0+), (+--0-)]
    """
    vectors = list(set(iterable))
    temp_vectors = vectors.copy()
    for i in irrelevant_components:
        for x in vectors:
            temp_vectors.append(x.set_to_plus([i]))
            temp_vectors.append(x.set_to_minus([i]))

        vectors = temp_vectors
        temp_vectors = vectors.copy()

    return vectors

def contraction(iterable: set[SignVector], indices: list[int]) -> set[SignVector]:
    r"""
    Return all sign vectors or vectors that are zero on given components.

    INPUT:

    - ``iterable`` -- an iterable of sign vectors or vectors

    - ``indices`` -- a list of indices.

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
    """
    return set(X.delete_components(indices) for X in iterable if not any(e in indices for e in X.support()))


def deletion(iterable: set[SignVector], indices: list[int]) -> set[SignVector]:
    r"""
    Remove given components from an iterable of sign vectors

    INPUT:

    - ``iterable`` -- an iterable of sign vectors

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
    """
    return set(X.delete_components(indices) for X in iterable)


def plot_sign_vectors(iterable: set[SignVector], vertex_size: int = 600, figsize: int = None, aspect_ratio=None):
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

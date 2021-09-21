#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.misc.flatten import flatten
from sage.modules.free_module_element import vector

from sign_vectors import SignVector, sign_vector, zero_sign_vector

# TODO: define closure
def closure(W, separate=False):
    r"""
    Computes the closure of a list of sign vectors.
    
    INPUT:
    
    - ``W`` -- a list of sign vectors
    
    - ``separate`` -- boolean (default: ``False``)
    
    OUTPUT:
    
    If ``separate`` is false, return the closure of ``W``. (default)
    
    If ``separate`` is true, separate the closure into lists, where each element
    has the same number of zero entries.
    
    .. NOTE::
    
       The sign vector $X$ is in the closure of a set of sign vectors $W$
       if there exists $Y \in W$ with $X \leq Y$.
    
    EXAMPLES:
    
    We consider a list consisting of only one sign vector::
    
        sage: from sign_vectors import sign_vector, closure
        sage: W = [sign_vector("+-0")]
        sage: W
        [(+-0)]
        sage: closure(W)
        [(000), (+00), (0-0), (+-0)]
    
    With the optional argument ``separate=True``, we can separate the resulting
    list into three lists.
    Each sign vector in such a list has the same number of zero entries::
    
        sage: closure(W, separate=True)
        [[(000)], [(+00), (0-0)], [(+-0)]]

    Now, we consider a list of three sign vectors::
    
        sage: W = [sign_vector("++-"), sign_vector("-00"), sign_vector("0--")]
        sage: W
        [(++-), (-00), (0--)]
        sage: closure(W)
        [(000), (+00), (-00), (0+0), (0-0), (00-), (++0), (+0-), (0+-), (0--), (++-)]
        sage: closure(W, separate=True)
        [[(000)],
         [(+00), (-00), (0+0), (0-0), (00-)],
         [(++0), (+0-), (0+-), (0--)],
         [(++-)]]
    
    TESTS::
        
        sage: closure([])
        []
    """
    if not W:
        return []
    n = W[0].length()
    F = [[zero_sign_vector(n)]]
    F_new = []
    for i in range(n):
        X = zero_sign_vector(n)
        X[i] = 1
        for Z in W:
            if X <= Z:
                F_new.append(X)
                break
        Y = zero_sign_vector(n)
        Y[i] = -1
        for Z in W:
            if Y <= Z:
                F_new.append(Y)
                break
    F.append(F_new)
    for i in range(1,n+1):
        F_new = []
        for X in F[1]: # X has always |supp(X)| = 1
            for Y in F[i]:
                if len(set(X.support() + Y.support())) == i+1: # TODO: utilize that the supports are sorted
                    Z = X.compose(Y)
                    if Z not in F_new: # is this necessary?
                        for V in W:
                            if Z <= V:
                                F_new.append(Z)
                                break
        if F_new == []:
            break
        else:
            F.append(F_new)
            if len(F_new) == 1:
                break
    if separate:
        return F
    else:
        return flatten(F)


def subvector(F, R):
    r"""
    Returns a function that returns a sign vector or vector consisting of entries not in ``R``.
    Used by ``contraction`` and ``deletion``.
    """
    if F == []:
        raise ValueError('List is empty.')
    n = len(list(F[0]))
    S = [e for e in range(n) if e not in R] # S = E\R

    if isinstance(F[0], SignVector):
        def vec(v):
            return sign_vector(v.list_from_positions(S))
    else:
        def vec(v):
            return vector(v.list_from_positions(S))
    return vec


def contraction(F, R, keep_components=False):
    r"""
    Returns all sign vectors that are zero on ``R``. Also works for real vectors.
    
    INPUT:
    
    - ``F`` -- a list of sign vectors, a list of vectors, or a matrix.
    
    - ``R`` -- a list of indices.
    
    - ``keep_components`` -- a boolean (default: ``False``).
    
    OUTPUT:
    
    - If ``keep_components`` is false, remove entries in ``R``. (default)
    
    - If ``keep_components`` is true, keep entries in ``R``.
    
    EXAMPLES:
    
        sage: from sign_vectors import sign_vector, contraction
        sage: W = [sign_vector("++0"), sign_vector("-00"), sign_vector("00+")]
        sage: W
        [(++0), (-00), (00+)]
    
    Only the third sign vector has a zero at the component with index ``0``.
    Removing this component leads to the following result::
    
        sage: contraction(W, [0])
        [(0+)]
        sage: contraction(W, [1])
        [(-0), (0+)]
        sage: contraction(W, [2])
        [(++), (-0)]
    
    The second sign vector has zeros at positions ``1`` and ``2``::
    
        sage: contraction(W, [1, 2])
        [(-)]
    
    We take the examples from before. With ``keep_components=True``, we keep the
    zero components of the appropriate sign vectors::
    
        sage: contraction(W, [0], keep_components=True)
        [(00+)]
        sage: contraction(W, [1], keep_components=True)
        [(-00), (00+)]
        sage: contraction(W, [2], keep_components=True)
        [(++0), (-00)]
        sage: contraction(W, [1, 2], keep_components=True)
        [(-00)]
    """
    if F == []:
        return F

    if keep_components:
        def vec(v):
            return v
    else:
        vec = subvector(F, R)
    
    L = []
    for X in F:
        val = True
        for e in X.support():
            if e in R:
                val = False
                break
        if val:
            L.append(vec(X))
    return L


def deletion(F, R):
    r"""
    Removes the components corresponding to ``R`` from a list of sign vectors.
    Also works for real vectors.
    
    INPUT:
    
    - ``F`` -- a list of sign vectors, a list of real vectors, or a matrix.
    
    - ``R`` -- a list of indices.
    
    EXAMPLES:
    
        sage: from sign_vectors import sign_vector, deletion                            
        sage: W = [sign_vector("+00"), sign_vector("++0"), sign_vector("00-")]          
        sage: W                                                                         
        [(+00), (++0), (00-)]
        sage: deletion(W, [0])                                                          
        [(00), (+0), (0-)]
    
    Duplicate sign vectors are removed if they would occur::
    
        sage: deletion(W, [1])                                                          
        [(+0), (0-)]
        sage: deletion(W, [1, 2])                                                       
        [(+), (0)]
    """
    if F == []:
        return F
    
    vec = subvector(F, R)
    L = []
    for X in F:
        X_R = vec(X)
        if X_R not in L: # naive, might be inefficient 
            L.append(X_R)
    return L


def loops(W):
    r"""
    Computes the list of loops of a given list of sign vectors (or real vectors)
    ``W``.
    
    ..NOTE::
    
        A loop is a component where every sign vector of ``W`` is zero.
    
    EXAMPLES::
    
        sage: from sign_vectors.utility import loops
        sage: from sign_vectors import sign_vector
        sage: loops([sign_vector([0,1,0]), sign_vector([-1,0,0])])
        [2]
        sage: loops([sign_vector([1,0,0]), sign_vector([-1,0,0])])
        [1, 2]
    
    Also works for real vectors::
    
        sage: loops([vector([5,0,0]), vector([2,0,-3])])
        [1]
    """
    if not W:
        raise ValueError('List is empty.')
    n = W[0].length()
    
    L = []
    for e in range(n):
        val = True
        for X in W:
            if X[e] != 0:
                val = False
                break
        if val == True:
            L.append(e)
    
    return L


def is_parallel(W, e, f, return_ratio=False):
    r"""
    Determines whether two elements ``e, f`` are parallel for each vector of ``W``.
    This also works for a set of sign vectors.
    
    INPUT:
    
    - ``W`` -- a list of vectors or sign vectors of length ``n``
    
    - ``e`` -- an integer with ``0 <= e <= n-1``

    - ``f`` -- an integer with ``0 <= e <= n-1``

    - ``return_ratio`` -- a boolean (default: False)

    OUTPUT:
    
    Returns a boolean.
    If ``return_ratio`` is true, a list consisting of the boolean and the ratio will be returned instead. 

    .. NOTE::
    
        The elements ``e`` and ``f`` are parallel if there exists a ratio ``d`` such that
        ``v[e] = d v[f]`` for each ``v`` in ``W``.
    
    EXAMPLES::
    
        sage: from sign_vectors.utility import is_parallel
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: is_parallel(L, 0, 1)
        False
        sage: is_parallel(L, 1, 2)
        False
        sage: is_parallel(L, 1, 3)
        True

    Now, we consider some real vectors::
    
        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: is_parallel(L, 0, 2)
        True
        sage: is_parallel(L, 0, 1)
        False
        sage: is_parallel(L, 1, 3)
        False
        sage: is_parallel(L, 4, 5)
        True
    
    We can also return the ratio of the two components::

        sage: is_parallel(L, 0, 2, return_ratio=True)
        [True, 1/2]
        sage: is_parallel(L, 2, 0, return_ratio=True)
        [True, 2]
        sage: is_parallel(L, 0, 1, return_ratio=True)                                   
        [False, 0]

    Also works for matrices::
    
        sage: M = matrix([[0,0,1,-1,0],[1,0,0,0,1],[1,1,1,1,1]])
        sage: is_parallel(M,0,4)
        True
        sage: is_parallel(M,0,1)
        False
    """
    d = 0 # will later be set to the ratio of X[e] and X[f]
    
    if return_ratio:
        ret = [False, 0]
    else:
        ret = False
    
    for X in W:
        if d == 0:
            if X[f] == 0:
                if X[e] != 0:
                    return ret
            elif X[e] == 0:
                return ret
            else: # determine ratio
                d = X[e]/X[f]
        else:
            if X[e] != d * X[f]:
                return ret
    if return_ratio:
        return [True, d]
    else:
        return True


def parallel_classes(W, positive_only=False):
    r"""
    Computes the parallel classes of a given set of vectors ``W``.
    This also works for a set of sign vectors.

    INPUT:

    - ``W`` -- a list of vectors or sign vectors of length ``n``
    
    - ``positive_only`` -- a boolean (default: False)
    
    OUTPUT:
    
    Returns a partition of ``[0, ..., n-1]`` into parallel classes.
    
    If ``positive_only`` is true, returns a partition of ``[0, ..., n-1]`` into positive parallel classes,
    that is, the ratios of the corresponding classes are non-negative.
    
    .. NOTE::
    
        The elements ``e`` and ``f`` are parallel if there exists a ratio ``d`` such that
        ``v[e] = d v[f]`` for each ``v`` in ``W``.
    
    EXAMPLES::
    
        sage: from sign_vectors.utility import parallel_classes
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: parallel_classes(L)
        [[0], [1, 3], [2]]
        sage: parallel_classes(L, positive_only=True)
        [[0], [1], [2], [3]]
    
    Now, we compute the parallel classes of a list of real vectors::
    
        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: parallel_classes(L)
        [[0, 2], [1], [3], [4, 5]]
    
    Let us compute the parallel classes of the rows of a matrix::
    
        sage: M = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: parallel_classes(M)
        [[0, 4], [1], [2, 3]]
        sage: parallel_classes(M, positive_only=True)
        [[0, 4], [1], [2], [3]]
    """
    if not W:
        raise ValueError('List is empty.')
    L = []
    k = W[0].length()
    toCheck = list(range(k))
    
    if positive_only:
        def is_par(W, e, f):
            val = is_parallel(W, e, f, return_ratio=True)
            if val[0] == False:
                return False
            elif val[1] < 0:
                return False
            else:
                return True
    else:
        def is_par(W, e, f):
            return is_parallel(W, e, f)
    
    while len(toCheck) > 0:
        e = toCheck.pop(0)
        l = [e]
        # `toCheck` might change in the for loop. -> toCheck[:]
        for f in toCheck[:]: # find parallel class `l` of ``e``
            if is_par(W, e, f):
                l.append(f)
                toCheck.remove(f)
        L.append(l)
    return L


def positive_parallel_classes(W):
    r"""
    Computes the positive parallel classes of a given set of vectors ``W``.
    This also works for a set of sign vectors.

    .. seealso::
    
        :func:`<sign_vectors.utility.parallel_classes>`

    EXAMPLES::
    
        sage: from sign_vectors.utility import positive_parallel_classes
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("--0+"), sign_vector("00+0")]
        sage: L
        [(++0-), (--0+), (00+0)]
        sage: positive_parallel_classes(L)
        [[0, 1], [2], [3]]
    
    Now, we compute the positive parallel classes of a list of real vectors::
    
        sage: L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        sage: L
        [(1, 1, 2, 3, 0, 0), (-2, 1, -4, 3, 3, -17), (0, 1, 0, 1, 0, 0)]
        sage: positive_parallel_classes(L)
        [[0, 2], [1], [3], [4], [5]]
    
    Let us compute the positive parallel classes of the rows of a matrix::
    
        sage: M = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]])
        sage: M
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: positive_parallel_classes(M)
        [[0, 4], [1], [2], [3]]
    """
    return parallel_classes(W, positive_only=True)


def classes_same_support(W):
    r"""
    Computes the classes with same support of a given list of sign vectors.
    Also works for a list of real vectors.

    INPUT:

    - ``W`` -- a list of vectors of length ``n``.
    
    EXAMPLES::
    
        sage: from sign_vectors.utility import classes_same_support
        sage: from sign_vectors import sign_vector
        sage: L = [sign_vector("++0-"), sign_vector("+-0+"), sign_vector("-0+0")]
        sage: L
        [(++0-), (+-0+), (-0+0)]
        sage: classes_same_support(L)
        [[(++0-), (+-0+)], [(-0+0)]]
        sage: classes_same_support([vector([1,1,0,0]), vector([2,-3,0,0]), vector([0,1,0,0])])
        [[(1, 1, 0, 0), (2, -3, 0, 0)], [(0, 1, 0, 0)]]
    """
    L = []
    Lc = [] # checked supports
    for X in W:
        s = X.support()
        if s not in Lc: # support of s has not been checked yet
            L.append([X]) # append new list with X
            Lc.append(s)
        else: # class of X already exists
            # find respective class of X
            for Li in L:
                if s == Li[0].support():
                    Li.append(X)
                    break
    return L

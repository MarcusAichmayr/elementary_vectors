#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

def loops(W):
    r"""
    Computes the list of loops of a given list of sign vectors ``W``.
    
    ..NOTE::
    
        A loop is a component where every sign vector of ``W`` is zero.
    """
    assert W, 'List is empty.'
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
        sage: W = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]]); W
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: parallel_classes(W)
        [[0, 4], [1], [2, 3]]
        sage: parallel_classes(W,positive_only=True)
        [[0, 4], [1], [2], [3]]
        """
    assert W, 'List is empty.'
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
        sage: W = matrix([[0,0,1,-2,0],[1,0,0,0,1],[1,1,-3,6,1]]); W
        [ 0  0  1 -2  0]
        [ 1  0  0  0  1]
        [ 1  1 -3  6  1]
        sage: positive_parallel_classes(W)
        [[0, 4], [1], [2], [3]]
    """
    return parallel_classes(W, positive_only=True)


def classes_same_support(W):
    r"""
    Computes the classes with same support of a given list of sign vectors.

    INPUT:

    - ``W`` -- a list of vectors of length ``n``.
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

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


def is_parallel(W,e,f):
    r"""
    Determines whether ``e`` is parallel to ``f`` in ``W``.
    
    .. NOTE::
    
        Works for sign vectors and for real vectors.
    """
    d = 0 # will later be set to the ratio of X[e] and X[f]
    for X in W:
        if d == 0:
            if X[f] == 0:
                if X[e] != 0:
                    return False
            elif X[e] == 0:
                return False
            else: # determine ratio
                d = X[e]/X[f]
        else:
            if X[e] != d * X[f]:
                return False
    return True


def parallel_classes(W):
    r"""
    Computes the parallel classes of a given set of sign vectors ``W``.

    INPUT:

    - ``W`` -- a list of sign vectors or vectors of length ``n``.
    
    OUTPUT:
    
    - a partition of ``[0, ..., n-1]`` into parallel classes.
    """
    assert W, 'List is empty.'
    L = []
    k = W[0].length()
    toCheck = list(range(k))
    
    while len(toCheck) > 0:
        e = toCheck.pop(0)
        l = [e]
        # `toCheck` might change in the for loop. -> toCheck[:]
        for f in toCheck[:]: # find parallel class `l` of ``e``
            if is_parallel(W, e, f):
                l.append(f)
                toCheck.remove(f)
        L.append(l)
    return L


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

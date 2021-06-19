#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sign_vectors import zero_sign_vector

from sage.misc.flatten import flatten

# Todo: define closure
def closure(W, separate=False):
    r"""
    Computes the closure of a list of sign vectors.
    
    INPUT:
    
    - ``W`` -- a list of sign vectors
    
    - ``separate`` -- boolean (default: ``False``)
    
    OUTPUT:
    
    If ``separate`` is false (default), return the closure of ``W``.
    
    If ``separate`` is true, separate the closure into lists, where each element
    has the same number of zero entries.
    """
    assert W, 'W is empty.'
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
                if len(set(X.support() + Y.support())) == i+1: # Todo: utilize that the supports are sorted
                    Z = X.compose(Y)
                    if Z not in F_new: # notwendig?
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

#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sign_vectors import sign_vector
from sign_vectors import zero_sign_vector

from sign_vectors.utility import classes_same_support, parallel_classes

def lower_faces(W):
    r"""
    Computes a list of lower faces.
    
    INPUT:
    
    - ``W`` -- a list of all covectors with same rank ``r`` of an oriented matroid.
    
    OUTPUT:
    
    Returns a list of covectors of rank ``r-1`` of the oriented matroid.
    """
    assert W, 'List is empty.'
    n = W[0].length()
    L = classes_same_support(W)
    W_ = []
    for Wj in L:
        PC = parallel_classes(Wj)
        for X in Wj:
            for D in PC:
                for i in D:
                    if X[i] != 0: # hence X_D != 0
                        if X.reverse_signs_in(D) in Wj:
                            Z = sign_vector([0 if i in D else X[i] for i in range(n)])
                            if Z not in W_: # Is this useless? Is Z in W_ for Z != 0 possible? - Yes, B = matrix(3,5,[1,2,3,5,7,1,4,2,5,3,6,2,1,2,3])
                                W_.append(Z)
                        break
    return W_


def face_enumeration(W):
    r"""
    Computes all covectors with less rank than the given list of covectors.
    
    INPUT:
    
    - ``W`` -- a list of all covectors of same rank ``r`` of an oriented matroid.
    
    OUTPUT:
    
    Returns a list of lists. Every list consists of all covectors of the same rank
    smaller than or equal to ``r`` of the oriented matroid.
    """
    assert W, 'List is empty.'
    L = [W]
    n = W[0].length()
    i = 0
    while L[0] != [zero_sign_vector(n)] and L[0] != []:
        L.insert(0, lower_faces(L[0]))
    return L

#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.combination import Combinations
from sage.modules.free_module_element import vector, zero_vector
from sage.structure.element import get_coercion_model
from sage.functions.other import binomial
from sage.arith.misc import gcd
from .utility import reduce_by_support

def elementary_vectors(M, kernel=True, reduce=True):
    r"""
    Computes elementary vectors of a subspace determined by the matrix ``M``.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    - ``kernel`` -- a boolean (default: ``True``)
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    OUTPUT:
    
    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``M``. (default)

    - If ``kernel`` is false, returns a list of elementary vectors lying in
      the row space of ``M``.
    
    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
      
    EXAMPLES::
    
        sage: from elementary_vectors import elementary_vectors
        sage: M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors(M)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, kernel=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, kernel=False)
        [(0, 1, 2, 0, 0), (0, 1, 0, 2, 0), (1, 0, 0, 0, 1), (0, 0, 1, -1, 0)]
    
    By default, the output is reduced, that is, each returned elementary vector
    has distinct support and common factors are canceled.
    
        sage: elementary_vectors(M, reduce=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, reduce=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]
    """
    if kernel:
        try:
            ind = M.pivot_rows() # would not work for polynomial matrices
            M = M.matrix_from_rows(ind)
            full_rank = True
        except:
            full_rank = False # The matrix might not have full rank. In this case, we will obtain ``[]``.
    else:
        M = M.right_kernel_matrix() # not implemented for polynomial matrices

    n = M.ncols()
    r = M.nrows()
    m = M.minors(r)

    L = elementary_vectors_from_minors(m, [r, n], reduce=reduce, ring=M.base_ring())
    if L == [] and not full_rank:
        print('WARNING: Could not determine rank of matrix. Result might be wrong.')
    return L


def elementary_vectors_from_minors(m, dim, reduce=True, ring=None):
    r"""
    Computes elementary vectors determined by given maximal minors of a matrix.
    
    INPUT:
    
    - ``m`` -- a list of maximal minors of a matrix
    
    - ``dim`` -- the dimensions of the matrix corresponding to ``m``
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    - ``ring`` -- parent of the entries of the matrix (despite the name,
                  this is not a priori required to be a ring).
                  By default, determine this from the given entries.
    
    OUTPUT:
    
    Uses the maximal minors ``m`` to compute the elementary vectors of the
    corresponding matrix.

    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
    
    EXAMPLES::
    
        sage: from elementary_vectors import elementary_vectors_from_minors
        sage: M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: m = M.minors(3); m
        [2, -2, 0, -4, 0, 0, 0, 2, -2, -4]
        sage: elementary_vectors_from_minors(m,[3,5])
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_minors(m, M.dimensions())
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_minors(m, M.dimensions(), reduce=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]
        sage: elementary_vectors_from_minors(m, M.dimensions(), ring=QQ)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
    """
    L = []
    r, n = dim
    if n <= r:
        return L
    assert binomial(n,r) == len(m), 'Dimensions do not fit.'
    if m == []:
        return L
    C = Combinations(n,r+1)
    dets = Combinations(n,r)
    if ring == None:
        # see SageMath/local/lib/python3.8/site-packages/sage/matrix/args.pyx
        # find common parent of minors
        ring = get_coercion_model().common_parent(*m)
    S = []
    for I in C: # construct vectors
        v = zero_vector(ring, n)
        for k in I:
            pos = I.index(k)
            Jk = list(I)
            Jk.pop(pos)
            v[k] = (-1)**pos * m[dets.rank(Jk)]
        if v != zero_vector(n):
            L.append(v)

    if reduce and L != []:
        out = []
        for v in reduce_by_support(L):
            out.append(vector(ring, v/gcd(v)))
        return out
    else:
        return L

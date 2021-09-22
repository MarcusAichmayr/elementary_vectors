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
from sign_vectors import sign_vector
from sage.symbolic.ring import SR
from .reductions import reduce_vectors, reduce_vector_using_equalities
from .utility import sign_determined
from sage.symbolic.assumptions import assume, forget, assumptions
import warnings

def elementary_vectors(data, dim=None, kernel=True, reduce=True, return_minors=False, ring=None):
    r"""
    Computes elementary vectors of a subspace determined by the matrix or from
    a list of maximal minors.
    
    INPUT:
    
    - ``data`` -- a matrix or a list of maximal minors
    
    - ``dim`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.
    
    - ``kernel`` -- a boolean (default: ``True``)
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    - ``return_minors`` -- a boolean (default: ``False``)
    
    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``data``.
    
    OUTPUT:
    
    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``data``. (default)
      Otherwise, returns a list of elementary vectors lying in
      the row space of ``data``.
      This argument is ignored if ``data`` is not a matrix.
      
    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
    
    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.
    
    .. SEEALSO::
        
        :func: `elementary_vectors_from_matrix`
        :func: `elementary_vectors_from_minors`
    
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
        sage: elementary_vectors(M, return_minors=True)
        [[(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)], [1, -1, 0, -2, 0, 0, 0, 1, -1, -2]]

    By default, the output is reduced, that is, each returned elementary vector
    has distinct support and common factors are canceled::
    
        sage: elementary_vectors(M, reduce=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(M, reduce=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]

    Here, ``data`` is a list of maximal minors::
    
        sage: M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: m = M.minors(3); m
        [2, -2, 0, -4, 0, 0, 0, 2, -2, -4]
        sage: elementary_vectors(m, [3,5])
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(m, M.dimensions())
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors(m, M.dimensions(), reduce=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]
        sage: elementary_vectors(m, M.dimensions(), ring=QQ)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
    """
    if isinstance(data, list):
        if dim == None:
            raise TypeError("When computing elementary vectors from a list of " +
                            "maximal minors, the dimensions of the corresponding " +
                            "matrix are needed.")
        evs = elementary_vectors_from_minors(data, dim=dim, reduce=reduce, ring=ring)
        if return_minors:
            return [evs, data]
        else:
            return evs
    else:
        return elementary_vectors_from_matrix(data, kernel=kernel, reduce=reduce, return_minors=return_minors, ring=ring)


def elementary_vectors_from_matrix(M, kernel=True, reduce=True, return_minors=False, ring=None):
    r"""
    Computes elementary vectors of a subspace determined by the matrix ``M``.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    - ``kernel`` -- a boolean (default: ``True``)
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    - ``return_minors`` -- a boolean (default: ``False``)
    
    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``M``.
    
    OUTPUT:
    
    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``M``. (default)

    - If ``kernel`` is false, returns a list of elementary vectors lying in
      the row space of ``M``.
    
    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
    
    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.
        
    EXAMPLES::
    
        sage: from elementary_vectors.elementary_vectors import elementary_vectors_from_matrix
        sage: M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]]); M
        [ 0  0  1 -1  0]
        [ 2  0  0  0  2]
        [ 1  1  1  1  1]
        sage: elementary_vectors_from_matrix(M)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_matrix(M, kernel=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_matrix(M, kernel=False)
        [(0, 1, 2, 0, 0), (0, 1, 0, 2, 0), (1, 0, 0, 0, 1), (0, 0, 1, -1, 0)]
        sage: elementary_vectors_from_matrix(M, return_minors=True)
        [[(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)], [1, -1, 0, -2, 0, 0, 0, 1, -1, -2]]
    
    By default, the output is reduced, that is, each returned elementary vector
    has distinct support and common factors are canceled::
    
        sage: elementary_vectors_from_matrix(M, reduce=True)
        [(0, 2, -1, -1, 0), (1, 0, 0, 0, -1)]
        sage: elementary_vectors_from_matrix(M, reduce=False)
        [(0, 4, -2, -2, 0),
         (2, 0, 0, 0, -2),
         (-2, 0, 0, 0, 2),
         (-4, 0, 0, 0, 4),
         (0, -4, 2, 2, 0)]
    
    Example with ``return_minors=True`` # TODO
    """
    if kernel:
        try:
            ind = M.pivot_rows() # would not work for polynomial matrices
            M = M.matrix_from_rows(ind)
            full_rank = True
        except:
            full_rank = False # The matrix might not have full rank. In this case, we will obtain ``[]``.
    else:
        M = M.right_kernel_matrix() # not implemented for polynomial matrices over ZZ (QQ does work.)

    n = M.ncols()
    r = M.nrows()
    m = M.minors(r)
    g = gcd(m)
    if g == 0:
        if not full_rank:
            warnings.warn('Could not determine rank of matrix. Result might be wrong.')
            return []
    if reduce:
        m = [mi/g for mi in m]
    if ring == None:
        ring = M.base_ring()
    evs = elementary_vectors_from_minors(m, [r, n], reduce=reduce, ring=ring)
    if return_minors:
        return [evs, m]
    else:
        return evs


def elementary_vectors_from_minors(m, dim, reduce=True, ring=None):
    r"""
    Computes elementary vectors determined by given maximal minors of a matrix.
    
    INPUT:
    
    - ``m`` -- a list of maximal minors of a matrix
    
    - ``dim`` -- the dimensions of the matrix corresponding to ``m``
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``m``.
    
    OUTPUT:
    
    Uses the maximal minors ``m`` to compute the elementary vectors of the
    corresponding matrix.

    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
    
    EXAMPLES::
    
        sage: from elementary_vectors.elementary_vectors import elementary_vectors_from_minors
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
    evs = []
    r, n = dim
    if n <= r:
        return evs
    if binomial(n,r) != len(m):
        raise ValueError('Dimensions do not fit.')
    if m == []:
        return evs
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
            evs.append(v)

    if reduce and evs != []:
        return reduce_vectors(evs)
#        out = []
#        for v in reduce_by_support(evs):
#            out.append(vector(ring, v/gcd(v)))
#        return out
    else:
        return evs


def non_negative_vectors(L):
    r"""
    Returns non-negative vectors.
    
    INPUT:
    
    - ``L`` -- a list of vectors
    
    OUTPUT:
    
    Returns all vectors of ``L`` that are
    - non_negative in each component; or
    - negative in each component. Those will be multiplied by ``-1``; or
    - containing variables such that no opposing signs occur.
    
    EXAMPLES::
    
        sage: from elementary_vectors import non_negative_vectors
        sage: l = [vector([1,1,0,-1]), vector([0,0,0,0]), vector([1,0,0,1])]
        sage: l
        [(1, 1, 0, -1), (0, 0, 0, 0), (1, 0, 0, 1)]
        sage: non_negative_vectors(l)
        [(0, 0, 0, 0), (1, 0, 0, 1)]

    Now, we consider an example with a variable::
    
        sage: from elementary_vectors import elementary_vectors
        sage: var('a')
        a
        sage: A = matrix([[a,0,0,0,1],[0,1,0,0,1]])
        sage: evs = elementary_vectors(A)
        sage: evs
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (-1, -a, 0, 0, a)]
        sage: non_negative_vectors(evs)
        ...
        UserWarning: Cannot determine sign of symbolic expression, returning 0 instead.
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (1, a, 0, 0, -a)]
        sage: assume(a > 0)
        sage: non_negative_vectors(evs)
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]
    """
    out = []
    for v in L:
        if sign_vector(v) >= 0: # Use ``>=`` instead of ``>``, (0,0,x) -> (000) should be returned
            out.append(v)
        elif sign_vector(v) < 0:
            out.append(-v)
    return out


# TODO: should assumptions be an optional argument?
def positive_elementary_vectors(data, dim=None, kernel=True, reduce=True, return_minors=False, ring=None):
    r"""
    Computes positive elementary vectors.
    
    INPUT:
    
    - ``data`` -- a matrix or a list of maximal minors
    
    - ``dim`` -- Not needed if ``data`` is a matrix.
                 Otherwise, this is a list of dimensions of the matrix
                 corresponding to the list of maximal minors ``data``.
    
    - ``kernel`` -- a boolean (default: ``True``)
    
    - ``reduce`` -- a boolean (default: ``True``)
    
    - ``return_minors`` -- a boolean (default: ``False``)
    
    - ``ring`` -- Parent of the entries of the elementary vectors.
                  By default, determine this from ``data``.
    
    OUTPUT:
    
    The output is a list of pairs.
    Each pair consists of a list of assumptions and the corresponding positive elementary vectors of ``data``.
    
    - If ``kernel`` is true, returns a list of elementary vectors lying in
      the kernel of ``data``. (default)
      Otherwise, returns a list of elementary vectors lying in
      the row space of ``data``.
      This argument is ignored if ``data`` is not a matrix.
      
    - If ``reduce`` is true, returned elementary vectors have distinct support
      and common factors are canceled.
    
    - If ``return_minors`` is true, a list is returned where the first
      element is the list of elementary vectors and the second element is
      the list of maximal minors used to compute the elementary vectors.
      
    EXAMPLES::
    
        sage: from elementary_vectors import positive_elementary_vectors
        sage: A = matrix([[1,-1,0]])                                                    
        sage: positive_elementary_vectors(A)                                            
        [[[], [(1, 1, 0), (0, 0, 1)]]]
        sage: positive_elementary_vectors(A, return_minors=True)                        
        [[[], [(1, 1, 0), (0, 0, 1)], [1, -1, 0]]]
        sage: M = matrix([[0,0,1,-1,0],[2,0,0,0,2],[1,1,1,1,1]])
        sage: positive_elementary_vectors(M)
        [[[], []]]
        
        TODO: Do more examples also symbolic examples
    """
    kwargs = locals()
    kwargs["return_minors"] = True
    kwargs["reduce"] = False
    
    def rec(i, l, eq):
        if i < len(m):
            if not sign_determined(m[i]):
                a = SR(m[i])
                try:
                    expr = a > 0
                    assume(expr)
                    rec(i+1, l + [expr], eq)
                    forget(expr)
                except ValueError:
                    pass
                try:
                    expr = a < 0
                    assume(expr)
                    rec(i+1, l + [expr], eq)
                    forget(expr)
                except ValueError:
                    pass
                try:
                    expr = a == 0
                    assume(expr)
                    rec(i+1, l + [expr], eq + [expr, -a == 0]) # a.substitute(-a == 0) results in a instead of 0.
                    forget(expr)
                except ValueError:
                    pass
            else:
                rec(i+1,l, eq)
        else:
            m_eq = reduce_vector_using_equalities(vector(m), eq=eq)
            if m_eq == 0: # minors are all zero
                warnings.warn('For ' + str(l) + ' all maximal minors are zero. There might be missing positive elementary vectors.')
#                 if isinstance(data, list):
#                     out.append([l, 'maximal minors are zero, could not compute elementary vectors for this branch'])
#                     return
#                 else:
#                     print('Todo')
#                     # substitute in data
#                     # compute positive_elementary_vectors from resulting matrix, eventuell auch l übergeben.
#                     # append each pair of result to out, also append assumptions l to first arguments
#                     pass
            # Do not overwrite ``evs`` here! It might be used in another branch.
            if reduce:
                L = reduce_vectors(evs, eq=eq) # TODO: this might result in an empty list. Instead, compute elementary vectors again if data is matrix
            else:
                L = evs
            if return_minors:
                out.append([l, non_negative_vectors(L), list(m_eq)])
            else:
                out.append([l, non_negative_vectors(L)])
            return
    
    evs, m = elementary_vectors(**kwargs)
    if evs == []:
        return []
    out = []
    rec(0, [], [])
    return out

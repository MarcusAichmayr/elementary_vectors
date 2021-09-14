#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.functions.generalized import sign
from sage.modules.free_module_element import zero_vector
from sage.misc.prandom import randint
from sage.rings.integer_ring import ZZ
from sage.symbolic.ring import SR
import warnings

class SignVector(SageObject):
    r"""A sign vector is an element of ``{-,+,0}^n``."""
    def __init__(self, l):
        r"""Creates a sign vector from a list ``l``."""
        try:
            self.__sv = vector(ZZ, [sign(x) for x in l])
        except:
            def sign_sym(a):
                r"""Returns appropriate sign of symbolic expression. Prints warning and returns ``0`` if sign cannot be computed."""
                if SR(a) > 0:
                    return 1
                elif SR(a) < 0:
                    return -1
                elif SR(a) == 0:
                    return 0
                else:
                    warnings.warn('Cannot determine sign of symbolic expression, returning 0 instead.')
                    return 0
            self.__sv = vector(ZZ, [sign_sym(x) for x in l])

    def _repr_(self):
        r"""Represents a sign vector by a string containing '-', '+' and '0'."""
        return '(' + ''.join([('+' if x > 0 else ('-' if x < 0 else '0')) for x in self.__sv]) + ')'
    
    def __hash__(self):
        return hash(repr(self))

    def length(self):
        r"""
        Returns the length of the sign vector.
        
        Examples::
        
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: X.length()
            3
        """
        return self.__sv.length()
    
    def __len__(self):
        r"""
        Returns the length of the sign vector.
        
        Examples::
        
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: len(X)
            3
        """
        return self.length()
    
    def compose(left, right):
        r"""
        Returns the composition of two sign vectors.
        
        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:
        
        Composition of this sign vector with ``other``.
        
        .. NOTE::
        
            Alternatively, the operator ``&`` can be used.
        
        EXAMPLES::
        
            sage: X = sign_vector('+00'); X
            (+00)
            sage: Y = sign_vector([-1,-1,0]); Y
            (--0)
            sage: X.compose(Y)
            (+-0)
            sage: Y.compose(X)
            (--0)
            sage: Y & X
            (--0)

        """
        if left.length() != right.length():
            raise ValueError('Sign vectors have different length.')
        return sign_vector([right[i] if left[i] == 0 else left[i] for i in range(left.length())])
    
    def __and__(left, right):
        r"""
        Returns the composition of two sign vectors.
        
        .. SEEALSO::
        
            :meth: `compose`
        """
        return left.compose(right)

    def __mul__(self, other):
        r"""Multiplication with a scalar."""
        return sign_vector(self.__sv.__mul__(sign(other)))
    
    def __rmul__(self, other):
        r"""Right multiplication with a scalar."""
        return self*other
    
    def __neg__(self):
        r"""
        Returns the sign vectors multiplied by ``-1``.
        
        Examples::
        
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: -X
            (0-+)
        """
        return sign_vector(-self.__sv)
        
    def __getitem__(self, e):
        r"""Returns the element at position ``e`` of the sign vector."""
        if isinstance(e, slice):
            return sign_vector(self.__sv[e])
        else:
            return self.__sv[e]
    
    def __setitem__(self, e, a):
        r"""Sets the element at position ``e`` to ``sign(a)``."""
        self.__sv[e] = sign(a)

    def support(self):
        r"""
        Returns a list of indices where the sign vector is non-zero.
        
        EXAMPLES::
        
            sage: X = sign_vector([-1,0,1,-1,0]); X
            (-0+-0)
            sage: X.support()
            [0, 2, 3]
        """
        return self.__sv.support()
    
    def __s_support(self, s):
        r"""Returns a list of entries where the sign vector equals ``s``."""
        return [e for e in range(self.length()) if self[e] == s]
    
    def zero_support(self):
        r"""
        Returns a list of indices where the sign vector is zero.
        
        EXAMPLES::
        
            sage: X = sign_vector([-1,0,1,-1,0]); X
            (-0+-0)
            sage: X.zero_support()
            [1, 4]
        """
        return self.__s_support(0)
    
    def positive_support(self):
        r"""
        Returns a list of indices where the sign vector is positive.
        
        EXAMPLES::
        
            sage: X = sign_vector([-1,0,1,-1,0]); X
            (-0+-0)
            sage: X.positive_support()
            [2]
        """
        return self.__s_support(1)
    
    def negative_support(self):
        r"""
        Returns a list of indices where the sign vector is negative.
        
        EXAMPLES::
        
            sage: X = sign_vector([-1,0,1,-1,0]); X
            (-0+-0)
            sage: X.negative_support()
            [0, 3]
        """
        return self.__s_support(-1)
    
    def list_from_positions(self, S):
        r"""Returns a list of the entries in the list ``S``."""
        return self.__sv.list_from_positions(S)
    
    def is_vector(self):
        r"""Returns ``False`` since sign vectors are not vectors."""
        return False

    def separating_elements(self, other):
        r"""
        Computes the list of separating elements of two sign vectors.
        
        INPUT:
        
        - ``other`` -- sign vector
        
        OUTPUT:
        List of elements ``e`` such that ``self[e] == -other[e] != 0``.
        
        EXAMPLES::
        
            sage: X = sign_vector('++00-'); X
            (++00-)
            sage: Y = sign_vector([1,-2,1,2,5]); Y
            (+-+++)
            sage: X.separating_elements(Y)
            [1, 4]
        """
        if self.length() != other.length():
            raise ValueError('Sign vectors have different length.')
        return [e for e in self.support() if self[e] == -other[e]]
    
    def reverse_signs_in(self, S):
        r"""
        Reverses sign of given entries.
        
        INPUT:
        
        - ``S`` -- list of indices
        
        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``S`` are multiplied by ``-1``.
        """
        return sign_vector([-self[e] if e in S else self[e] for e in range(self.length())])

    def conforms(left, right):
        r"""
        Conformal relation of two sign vectors.
        
        EXAMPLES::
        
            sage: X = sign_vector([-1, 1, 0, 0, 1]); X
            (-+00+)
            sage: Y = sign_vector([-1, 1, 1, 0, 1]); Y
            (-++0+)
            sage: Z = sign_vector([-1, 1, 1, -1, 1]); Z
            (-++-+)
            sage: X.conforms(Y)
            True
            sage: X.conforms(X)
            True
            sage: Y.conforms(Z)
            True
            sage: Z.conforms(Y)
            False
            sage: X.conforms(Z)
            True
        """
        if left.length() != right.length():
            raise ValueError('Sign vectors have different length.')
        
        def lessthan(x,y):
            r"""Unary conformal relation."""
            if x == 0:
                return True
            elif x == y:
                return True
            else:
                return False
    
        for e in range(left.length()):
            if not lessthan(left[e], right[e]):
                return False
        return True
    
    def __eq__(self, other):
        r"""Returns whether this sign vector is equal to ``other``."""
        if isinstance(other, SignVector):
            return self.__sv == other.__sv
        else:
            return self.__sv == other
    
    def __le__(left, right):
        r"""Returns whether this sign vector is less or equal to ``right``."""
        if isinstance(right, SignVector):
            return left.conforms(right)
        elif right == 0:
            for a in left:
                if a > 0:
                    return False
            return True
        else:
            return left.__sv <= right # should this be that way?

    def __lt__(left, right):
        r"""Returns whether this sign vector is less than ``right``."""
        return left != right and left <= right
    
    def __ge__(left, right):
        r"""Returns whether this sign vector is greater or equal to ``right``."""
        if isinstance(right, SignVector):
            return right.conforms(left)
        elif right == 0:
            for a in left:
                if a < 0:
                    return False
            return True
        else:
            return left.__sv >= right
        
    def __gt__(left, right):
        r"""Returns whether this sign vector is greater than ``right``."""
        return left != right and left >= right
    
    def is_orthogonal_to(self, other):
        r"""
        Returns whether two sign vectors are orthogonal.
        
        INPUT:
        
        - ``other`` -- a sign vector.
        
        OUTPUT:
        
        - Returns ``True`` if the sign vectors are orthogonal and ``False`` otherwise.
        
        EXAMPLES::
        
            sage: X = sign_vector('-+00+'); X
            (-+00+)
            sage: X.is_orthogonal_to(X)
            False
            sage: X.is_orthogonal_to(sign_vector('++000'))
            True
            sage: X.is_orthogonal_to(sign_vector('++00+'))
            True
            sage: X.is_orthogonal_to(sign_vector('00++0'))
            True
        """
        if self.length() != other.length():
            raise ValueError('Sign vectors have different length.')

        if [e for e in self.support() if e in other.support()] == []:
            return True
        else:
            for e in self.support():
                if self[e]*other[e] > 0:
                    for f in self.support():
                        if self[f]*other[f] < 0:
                            return True
            return False

def sign_vector(v):
    r"""
    Create a sign vector from a list, vector or string.
    
    INPUT:
    
    - ``v`` -- different inputs are accepted:
    
        - an iterable (e.g. a list or vector) of real values.
          Variables can also occur.
    
        - a string consisting of '-', '+', '0'. Other characters are treated as '0'.
    
    OUTPUT:
    
    Returns a sign vector. If variables occur and the sign of the corresponding
    entries is not determined, prints a warning and inserts ``0`` instead.
    
    EXAMPLES::
    
        sage: from sign_vectors import sign_vector
        sage: sign_vector([5,0,-1,-2])
        (+0--)
        sage: v = vector([5,0,-1,-2])
        sage: sign_vector(v)
        (+0--)
      
    We can also use a string to compute a sign vector::
    
        sage: sign_vector('++-+-00-')
        (++-+-00-)]
        
    Variables are supported to some extent::
    
        sage: v = vector([1,x,-1])
        sage: V = sign_vector(v); V
        Warning: Cannot determine sign of symbolic expression, returning 0 instead.
        (+0-)
        sage: vector(V)
        (1, 0, -1)
        sage: assume(x > 0)
        sage: sign_vector(v)
        (++-)
    """
    if isinstance(v, str):
        return SignVector([1 if t == '+' else (-1 if t == '-' else 0) for t in v])
    else:
        return SignVector(list(v))

def zero_sign_vector(n):
    r"""
    Return the zero sign vector of length ``n``.
    
    Examples::
    
        sage: zero_sign_vector(5)
        (00000)
    """
    return sign_vector(zero_vector(n))

def random_sign_vector(n):
    r"""
    Return a random sign vector of length ``n``.#
    
    Examples::
    
        sage: random_sign_vector(5) # random
        (++-0-)
    """
    return sign_vector([randint(-1,1) for k in range(n)])

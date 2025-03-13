r"""
Sign vectors
============

First, we load the functions from the package::

    sage: from sign_vectors import *

There are several ways to define sign vectors::

    sage: sign_vector([1, 0, -1, 1])
    (+0-+)
    sage: x = vector([5, 2/5, -1, 0])
    sage: sign_vector(x)
    (++-0)
    sage: sign_vector('++-+-00-')
    (++-+-00-)

We can easily construct the zero sign vector of a given length::

    sage: zero_sign_vector(6)
    (000000)

It is also possible to generate some random sign vector::

    sage: random_sign_vector(7) # random
    (+-+00+-)

To reverse signs, we apply the operator ``-`` as usual::

    sage: X = sign_vector('+++000--')
    sage: X
    (+++000--)
    sage: -X
    (---000++)

There are different notions of support::

    sage: X.support()
    [0, 1, 2, 6, 7]
    sage: X.zero_support()
    [3, 4, 5]
    sage: X.positive_support()
    [0, 1, 2]
    sage: X.negative_support()
    [6, 7]

Next, we define two sign vectors::

    sage: X = sign_vector([-1, 0, 1, -1, 0])
    sage: X
    (-0+-0)
    sage: Y = sign_vector([0, 1, 0, 1, 0])
    sage: Y
    (0+0+0)

We can compose them::

    sage: X.compose(Y)
    (-++-0)
    sage: Y.compose(X)
    (-+++0)

One can also use the operator ``&`` to compose sign vectors::

    sage: X & Y
    (-++-0)
    sage: Y & X
    (-+++0)

The conformal relation is a partial order on a set of sign vectors::

    sage: X = sign_vector([-1, 1, 0, 0, 1])
    sage: Y = sign_vector([-1, 1, 1, 0, 1])
    sage: Z = sign_vector([-1, 1, 1, -1, 1])
    sage: X
    (-+00+)
    sage: Y
    (-++0+)
    sage: Z
    (-++-+)

We can apply it in the following way::

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
    sage: X < X
    False
    sage: X <= Y
    True
    sage: Y > Z
    False
    sage: Z >= X
    True

Similar as for real vectors, we define orthogonality for sign vectors.
First, we define some real vectors::

    sage: x = vector([0, 0, 1])
    sage: y = vector([1, -2, 0])
    sage: z = vector([2, 1, 2])

We compute some scalar products to investigate orthogonality of those vectors::

    sage: x.dot_product(y)
    0
    sage: y.dot_product(z)
    0
    sage: x.dot_product(z)
    2

Next, we define the corresponding sign vectors::

    sage: X = sign_vector(x)
    sage: X
    (00+)
    sage: Y = sign_vector(y)
    sage: Y
    (+-0)
    sage: Z = sign_vector(z)
    sage: Z
    (+++)

By definition, if two real vectors are orthogonal, then their corresponding
sign vectors are also orthogonal::

    sage: X.is_orthogonal_to(Y)
    True
    sage: Y.is_orthogonal_to(Z)
    True
    sage: X.is_orthogonal_to(Z)
    False
"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import warnings
from random import choices
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SR


length_error = ValueError("Elements have different length.")


def sign_symbolic(value) -> int:
    r"""
    Return sign of expression. Supports symbolic expressions.

    OUTPUT:
    If the sign cannot be determined, a warning is shown and ``0`` is returned.

    EXAMPLES::

        sage: from sign_vectors.sign_vectors import sign_symbolic
        sage: sign_symbolic(1)
        1
        sage: sign_symbolic(-2)
        -1
        sage: var('a')
        a
        sage: sign_symbolic(a)
        ...
        UserWarning: Cannot determine sign of symbolic expression, using ``0`` instead.
        0
        sage: assume(a > 0)
        sage: sign_symbolic(a)
        1
        sage: forget()
        sage: var('a, b, c, d')
        (a, b, c, d)
        sage: assume(a > 0, b > 0, c > 0, d > 0)
        sage: sign_symbolic((a + b)*(c + d) - b*c)
        1
    """
    if value == 0:
        return 0
    if value > 0:
        return 1
    if value < 0:
        return -1

    expr = SR(value).simplify_full()
    if expr == 0:
        return 0
    if expr > 0:
        return 1
    if expr < 0:
        return -1

    warnings.warn("Cannot determine sign of symbolic expression, using ``0`` instead.")
    return 0


class Sign(SageObject):
    r"""An element in ``{+, -, 0}``."""

    def __init__(self, value) -> None:
        if isinstance(value, Sign):
            self._positive = value._positive
            return

        value = sign_symbolic(value)
        if value > 0:
            self._positive = True
        elif value < 0:
            self._positive = False
        else:
            self._positive = None

    def _repr_(self) -> str:
        if self._positive is None:
            return "0"
        if self._positive:
            return "+"
        return "-"

    def __hash__(self):
        if self._positive is None:
            return 0
        if self._positive:
            return 2
        return 1

    def is_positive(self) -> bool:
        r"""Return whether this sign is ``+``."""
        return self._positive

    def is_negative(self) -> bool:
        r"""Return whether this sign is ``-``."""
        return self._positive is False

    def is_zero(self) -> bool:
        r"""Return whether this sign is ``0``."""
        return self._positive is None

    def compose(self, other):
        r"""
        Return the composition of two signs.

        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:

        Composition of this sign vector with ``other``.

        .. NOTE::

            Alternatively, the operator ``&`` can be used.

        EXAMPLES::

            sage: from sign_vectors.sign_vectors import Sign
            sage: Sign(1).compose(Sign(0))
            +
            sage: Sign(0).compose(Sign(-1))
            -
            sage: Sign(1).compose(Sign(-2))
            +
            sage: Sign(1) & Sign(2)
            +
            sage: Sign(-1) & Sign(2)
            -
            sage: Sign(0) & Sign(0)
            0
        """
        if self._positive is None:
            return other
        return self

    def __and__(self, other):
        r"""
        Return the composition of two signs.

        .. SEEALSO::

            :meth: `compose`
        """
        return self.compose(other)

    def __mul__(self, other):
        r"""Multiplication with a scalar."""
        if isinstance(other, Sign):
            if other.is_zero() or self.is_zero():
                return Sign(0)
            if self.is_positive():
                return other
            return -other
        return self * Sign(other)

    def __rmul__(self, other):
        r"""Right multiplication with a scalar."""
        return self * other

    def __neg__(self):
        r"""Multiplication with ``-1``."""
        if self.is_positive():
            return Sign(-1)
        if self.is_negative():
            return Sign(1)
        return Sign(0)

    def __truediv__(self, other) -> int:
        r"""
        Division of two signs.

        EXAMPLES::

            sage: from sign_vectors.sign_vectors import Sign
            sage: Sign(1) / Sign(-2)
            -1
            sage: Sign(0) / Sign(1)
            0
        """
        return self.to_integer() // other.to_integer()

    def __eq__(self, other) -> bool:
        r"""
        Return whether this sign is equal to ``other``.

        EXAMPLES::

            sage: from sign_vectors.sign_vectors import Sign
            sage: Sign(1) == Sign(2)
            True
            sage: Sign(1) == Sign(-2)
            False
            sage: Sign(-1) == Sign(0)
            False
        """
        if isinstance(other, Sign):
            return self._positive == other._positive
        if other == 0:
            return self.is_zero()
        return False

    def __le__(self, other) -> bool:
        r"""
        Return whether this sign is less or equal to ``other``.

        EXAMPLES::

            sage: from sign_vectors.sign_vectors import Sign
            sage: Sign(1) <= Sign(1)
            True
            sage: Sign(1) <= Sign(-1)
            False
            sage: Sign(-1) <= Sign(1)
            False
            sage: Sign(0) <= Sign(-1)
            True
            sage: Sign(0) <= Sign(0)
            True
        """
        if isinstance(other, Sign):
            if self.is_zero():
                return True
            return self == other
        if other == 0:
            return not self.is_positive()
        return False

    def __lt__(self, other) -> bool:
        r"""
        Return whether this sign is less than ``other``.

        .. SEEALSO::

            :meth: `compose`

        EXAMPLES::

            sage: from sign_vectors.sign_vectors import Sign
            sage: Sign(0) < Sign(1)
            True
            sage: Sign(0) < Sign(-1)
            True
            sage: Sign(1) < Sign(1)
            False

        We can also use ``<`` to compare a sign vector with ``0``::

            sage: Sign(-1) < 0 #TODO: Do we want this?
            True
            sage: Sign(0) < 0
            False
            sage: 0 < Sign(1)
            True
        """
        if isinstance(other, Sign):
            return self != other and self <= other
        if other == 0:
            return self.is_negative()
        raise TypeError("unsupported operand")

    def __ge__(self, other) -> bool:
        r"""
        Return whether this sign is greater or equal to ``other``.

        .. SEEALSO::

            :meth: `conforms`
        """
        if isinstance(other, Sign):
            return other.conforms(self)
        if other == 0:
            return not self.is_negative()
        raise TypeError("unsupported operand")

    def __gt__(self, other) -> bool:
        r"""
        Return whether this sign vector is greater than ``other``.

        .. SEEALSO::

            :meth: `conforms`
        """
        if isinstance(other, Sign):
            return self != other and self >= other
        if other == 0:
            return self.is_positive()
        raise TypeError("unsupported operand")

    def to_integer(self):
        r"""Return the related integer."""
        if self.is_positive():
            return 1
        if self.is_negative():
            return -1
        return 0


class SignVector(SageObject):
    r"""A sign vector."""

    __slots__ = ("_support", "_positive_support", "_length")

    def __init__(self, support: frozenset[int], psupport: frozenset[int], length: int) -> None:
        r"""
        Create a sign vector object.

        INPUT:

        - ``support`` -- a frozenset that represents the support of this sign vector.

        - ``psupport`` -- a frozenset that represents the positive support of this sign vector.

        - ``length`` -- the length of this sign vector.

        .. NOTE::

            The ``psupport`` should be a subset of ``support``.
            For efficiency, this is not checked when creating a sign vector object.

        .. SEEALSO::

            :func: `~sign_vector`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: s = frozenset([0, 1, 3])
            sage: p = frozenset([1, 3])
            sage: SignVector(s, p, 4)
            (-+0+)
        """
        self._support = support
        self._positive_support = psupport
        self._length = length

    def _repr_(self) -> str:
        return "(" + str(self) + ")"

    def __str__(self) -> str:
        return "".join(
            "+"
            if e in self._positive_support
            else ("-" if e in self._support else "0")
            for e in range(self.length())
        )

    def length(self) -> int:
        r"""
        Return the length of the sign vector.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: X.length()
            3
        """
        return self._length

    def __len__(self) -> int:
        r"""
        Return the length of this sign vector.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: len(X)
            3
        """
        return self._length

    def support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is nonzero.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 0, 1, -1, 0])
            sage: X
            (-0+-0)
            sage: X.support()
            [0, 2, 3]
        """
        return list(self._support)

    def zero_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is zero.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 0, 1, -1, 0])
            sage: X
            (-0+-0)
            sage: X.zero_support()
            [1, 4]
        """
        return [e for e in range(self.length()) if not e in self._support]

    def positive_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is positive.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 0, 1, -1, 0])
            sage: X
            (-0+-0)
            sage: X.positive_support()
            [2]
        """
        return list(self._positive_support)

    def negative_support(self) -> list[int]:
        r"""
        Return a list of indices where the sign vector is negative.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 0, 1, -1, 0])
            sage: X
            (-0+-0)
            sage: X.negative_support()
            [0, 3]
        """
        return list(self._negative_support())

    def _negative_support(self) -> frozenset[int]:
        r"""Return the set corresponding to the negative support."""
        return self._support.symmetric_difference(self._positive_support)

    def __getitem__(self, e):
        r"""
        Return the element at position ``e`` of the sign vector.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector("0++-")
            sage: X
            (0++-)
            sage: X[0]
            0
            sage: X[1]
            1
            sage: X[3]
            -1
            sage: X[-1]
            -1
            sage: X[1:3] # todo: not implemented
            (++)

        TESTS::

            sage: X[-2]
            1
            sage: X[100]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if isinstance(e, slice):
            raise NotImplementedError("TODO")
        if e >= self.length() or e < -self.length():
            raise IndexError("index out of range")
        if e < 0:
            e %= self.length()
        if e in self._support:
            return 1 if e in self._positive_support else -1
        return 0

    def list_from_positions(self, S) -> list[int]:
        r"""
        Return a list of components that are in the list of indices ``S``.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 1, 0, 0, 1])
            sage: X
            (-+00+)
            sage: X.list_from_positions([0, 1, 4])
            [-1, 1, 1]
        """
        return [self[e] for e in S]

    def compose(self, other):
        r"""
        Return the composition of two sign vectors.

        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:

        Composition of this sign vector with ``other``.

        .. NOTE::

            Alternatively, the operator ``&`` can be used.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('+00')
            sage: X
            (+00)
            sage: Y = sign_vector([-1, -1, 0])
            sage: Y
            (--0)
            sage: X.compose(Y)
            (+-0)
            sage: Y.compose(X)
            (--0)
            sage: Y & X
            (--0)
            sage: X = sign_vector('000+++---')
            sage: Y = sign_vector('0+-0+-0+-')
            sage: X.compose(Y)
            (0+-+++---)
        """
        if self.length() != other.length():
            raise length_error

        support = self._support.union(other._support)
        psupport = frozenset(
            e
            for e in support
            if (e in self._positive_support)
            or (e in other._positive_support and not e in self._support)
        )

        return SignVector(support, psupport, self.length())

    def compose_harmonious(self, other):
        r"""
        Return the composition of two harmonious sign vectors.

        INPUT:

        - ``other`` -- a sign vector

        OUTPUT:

        Composition of this sign vector with ``other``.

        .. NOTE::

            This method is more efficient than :meth:`compose`.
            However, it does not check whether the sign vectors are harmonious.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('+00')
            sage: X
            (+00)
            sage: Y = sign_vector([0, -1, 0])
            sage: Y
            (0-0)
            sage: X.compose_harmonious(Y)
            (+-0)
            sage: Y.compose_harmonious(X)
            (+-0)
        """
        if self.length() != other.length():
            raise length_error

        return SignVector(
            self._support.union(other._support),
            self._positive_support.union(other._positive_support),
            self.length(),
        )

    def __and__(self, other):
        r"""
        Return the composition of two sign vectors.

        .. SEEALSO::

            :meth: `compose`

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('+00')
            sage: X
            (+00)
            sage: Y = sign_vector([-1, -1, 0])
            sage: Y
            (--0)
            sage: X & Y
            (+-0)
            sage: Y & X
            (--0)
        """
        return self.compose(other)

    def __mul__(self, value):
        r"""
        Multiplication with a scalar.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 1, 0, 0, 1])
            sage: X
            (-+00+)
            sage: -1*X
            (+-00-)
            sage: 1*X
            (-+00+)
        """
        if value > 0:
            return +self
        if value < 0:
            return -self
        return zero_sign_vector(self.length())

    def __rmul__(self, value):
        r"""
        Right multiplication with a scalar.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 1, 0, 0, 1])
            sage: X
            (-+00+)
            sage: X*(-1)
            (+-00-)
            sage: X*1
            (-+00+)
        """
        return self * value

    def separating_elements(self, other) -> list[int]:
        r"""
        Compute the list of separating elements of two sign vectors.

        INPUT:

        - ``other`` -- sign vector

        OUTPUT:
        List of elements ``e`` such that ``self[e] == -other[e] != 0``.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('++00-')
            sage: X
            (++00-)
            sage: Y = sign_vector([1, -2, 1, 2, 5])
            sage: Y
            (+-+++)
            sage: X.separating_elements(Y)
            [1, 4]
        """
        if self.length() != other.length():
            raise length_error
        return [
            e
            for e in self._support.intersection(other._support)
            if (e in self._positive_support) ^ (e in other._positive_support)
        ]

    def is_harmonious(self, other) -> bool:
        r"""
        Check whether these two sign vectors are harmonious.

        INPUT:

        - ``other`` -- a sign vector or real vector

        OUTPUT:
        Returns true if there are no separating elements.
        Otherwise, false is returned.

        .. NOTE::

            Two sign vectors are harmonious if there is no component
            where one sign vector has ``+`` and the other has ``-``.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('++00-')
            sage: X
            (++00-)
            sage: Y = sign_vector([1, -2, 1, 2, 5])
            sage: Y
            (+-+++)
            sage: X.is_harmonious(Y)
            False
            sage: sign_vector('0+00').is_harmonious(sign_vector('-+0+'))
            True
            sage: v = vector([1, 2/3, 0, -1, -1])
            sage: X.is_harmonious(v)
            True
            sage: Y.is_harmonious(v)
            False
        """
        if self.length() != other.length():
            raise length_error
        if not isinstance(other, SignVector):
            other = sign_vector(other)

        for e in self._support:
            if e not in other._support:
                continue
            if (e in self._positive_support) ^ (e in other._positive_support):
                return False
        return True

    def disjoint_support(self, other) -> bool:
        r"""
        Return whether these two sign vectors have disjoint support.

        INPUT:

        - ``other`` -- a sign vector or real vector

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('++00')
            sage: Y = sign_vector('0+0-')
            sage: Z = sign_vector('00--')
            sage: X.disjoint_support(Y)
            False
            sage: Y.disjoint_support(X)
            False
            sage: X.disjoint_support(Z)
            True
            sage: Z.disjoint_support(X)
            True
            sage: Y.disjoint_support(Z)
            False
        """
        if self.length() != other.length():
            raise length_error
        if isinstance(other, SignVector):
            return self._support.isdisjoint(other._support)
        return self._support.isdisjoint(other.support())

    def reverse_signs_in(self, indices):
        r"""
        Reverses sign of given entries.

        INPUT:

        - ``indices`` -- list of indices

        OUTPUT:
        Returns a new sign vector of same length. Components with indices in
        ``indices`` are multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector([-1, 1, 1, 0, 1])
            sage: X
            (-++0+)
            sage: X.reverse_signs_in([0, 2, 3])
            (++-0+)
        """
        support = frozenset(self._support)
        psupport = set(self._positive_support)
        for e in indices:
            if e in support:
                if e in psupport:
                    psupport.remove(e)
                else:
                    psupport.add(e)
        return SignVector(support, frozenset(psupport), self.length())

    def conforms(self, other) -> bool:
        r"""
        Conformal relation of two sign vectors.

        .. NOTE::

            Alternatively, the operator ``<=`` can be used.
            Use ``>=``, ``<`` and ``>`` for the other relations.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
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
        if self.length() != other.length():
            raise length_error

        return all(self[e] == other[e] for e in self._support)

    def __eq__(self, other) -> bool:
        r"""
        Return whether this sign vector is equal to ``other``.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector("++0-")
            sage: X == X
            True
            sage: X == sign_vector("00++")
            False

        TESTS::

            sage: from sign_vectors import *
            sage: zero_sign_vector(3) == 0
            True
            sage: 0 == zero_sign_vector(3)
            True
        """
        if isinstance(other, SignVector):
            if self.length() != other.length():
                raise length_error
            return (
                self._support == other._support
                and self._positive_support == other._positive_support
            )
        if other == 0:
            return not self._support
        return False

    def __le__(self, other) -> bool:
        r"""
        Return whether this sign vector is less or equal to ``other``.

        .. SEEALSO::

            :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector([-1, 1, 0, 0, 1]); X
            (-+00+)
            sage: Y = sign_vector([-1, 1, 1, 0, 1]); Y
            (-++0+)
            sage: X <= Y
            True

        We can also use ``<=`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') <= 0
            True
            sage: sign_vector([1, 1, -1, 0]) <= 0
            False
            sage: 0 <= sign_vector([1, 1, 0, 0])
            True
            sage: zero_sign_vector(2) <= 0
            True
        """
        if isinstance(other, SignVector):
            return self.conforms(other)
        if other == 0:
            return all(a <= 0 for a in self)
        return False

    def __lt__(self, other) -> bool:
        r"""
        Return whether this sign vector is less than ``other``.

        .. SEEALSO::

            :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector([-1, 1, 0, 0, 1]); X
            (-+00+)
            sage: Y = sign_vector([-1, 1, 1, 0, 1]); Y
            (-++0+)
            sage: X < Y
            True

        We can also use ``<`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') < 0
            True
            sage: sign_vector([1, 1, -1, 0]) < 0
            False
            sage: 0 < sign_vector([1, 1, 0, 0])
            True
            sage: zero_sign_vector(2) < 0
            False
        """
        return self != other and self <= other

    def __ge__(self, other) -> bool:
        r"""
        Return whether this sign vector is greater or equal to ``other``.

        .. SEEALSO::

            :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector([-1, 1, 0, 0, 1]); X
            (-+00+)
            sage: Y = sign_vector([-1, 1, 1, 0, 1]); Y
            (-++0+)
            sage: Y >= X
            True

        We can also use ``>=`` to compare a sign vector with ``0``::

            sage: sign_vector('00--') >= 0
            False
            sage: sign_vector([1, 1, -1, 0]) >= 0
            False
            sage: sign_vector([1, 1, 0, 0]) >= 0
            True
            sage: zero_sign_vector(2) >= 0
            True
        """
        if isinstance(other, SignVector):
            return other.conforms(self)
        if other == 0:
            return all(a >= 0 for a in self)
        return False

    def __gt__(self, other) -> bool:
        r"""
        Return whether this sign vector is greater than ``other``.

        .. SEEALSO::

            :meth: `conforms`

        EXAMPLES::

            sage: from sign_vectors import *
            sage: X = sign_vector([-1, 1, 0, 0, 1]); X
            (-+00+)
            sage: Y = sign_vector([-1, 1, 1, 0, 1]); Y
            (-++0+)
            sage: Y > X
            True

        We can also use ``>`` to compare a sign vector with ``0``::

            sage: 0 > sign_vector('00--')
            True
            sage: sign_vector([1, 1, -1, 0]) > 0
            False
            sage: sign_vector([1, 1, 0, 0]) > 0
            True
            sage: zero_sign_vector(2) > 0
            False
        """
        return self != other and self >= other

    def is_orthogonal_to(self, other) -> bool:
        r"""
        Return whether two sign vectors are orthogonal.

        INPUT:

        - ``other`` -- a sign vector.

        OUTPUT:

        - Returns ``True`` if the sign vectors are orthogonal and ``False`` otherwise.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
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
            raise length_error
        if self._support.isdisjoint(other._support):
            return True
        have_positive_product = False
        have_negative_product = False
        for e in self._support:
            if self[e] * other[e] > 0:
                have_positive_product = True
            elif self[e] * other[e] < 0:
                have_negative_product = True
            if have_positive_product and have_negative_product:
                return True
        return False

    def __bool__(self) -> bool:
        return self != 0

    def __neg__(self):
        r"""
        Return the sign vectors multiplied by ``-1``.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('0+-'); X
            (0+-)
            sage: -X
            (0-+)
        """
        return SignVector(self._support, self._negative_support(), self.length())

    def __pos__(self):
        r"""
        Return the sign vectors multiplied by ``1`` which is a copy of this sign vector.

        EXAMPLES::

            sage: from sign_vectors import sign_vector
            sage: X = sign_vector('0+-')
            sage: X
            (0+-)
            sage: +X
            (0+-)
        """
        return SignVector(self._support, self._positive_support, self.length())

    def is_vector(self) -> bool:
        r"""Return ``False`` since sign vectors are not vectors."""
        return False

    def __hash__(self):
        r"""Return the hash value of this sign vector."""
        # TODO remove length
        return hash((self._length, self._support, self._positive_support))

    @staticmethod
    def from_str(s: str):
        r"""
        Creates a sign vector from a string.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.from_str("+-0+0")
            (+-0+0)
        """
        return SignVector(
            frozenset(pos for pos, t in enumerate(s) if t in "+-"),
            frozenset(pos for pos, t in enumerate(s) if t == "+"),
            len(s),
        )

    @staticmethod
    def from_support(support: list, psupport: list, length: int):
        r"""
        Return a sign vector that is given by lists representing support and positive support.

        INPUT:

        - ``support`` -- a list

        - ``psupport`` -- a list

        - ``length`` -- a nonnegative integer

        OUTPUT:
        a sign vector

        .. NOTE::

            The list ``psupport`` should be a sublist of ``support``.
            For efficiency, this is not checked.

        EXAMPLES::

            sage: from sign_vectors import *
            sage: SignVector.from_support([1, 2, 4], [1, 4], 6)
            (0+-0+0)
        """
        return SignVector(frozenset(support), frozenset(psupport), length)


def sign_vector(iterable):
    r"""
    Create a sign vector from a list, vector or string.

    INPUT:

    - ``iterable`` -- different inputs are accepted:

        - an iterable (e.g. a list or vector) of real values.
          Variables can also occur.

        - a string consisting of ``"-"``, ``"+"``, ``"0"``. Other characters are treated as ``"0"``.

    OUTPUT:

    Returns a sign vector. If variables occur and the signs of the corresponding
    entries cannot be determined, prints a warning and inserts ``"0"`` instead.

    EXAMPLES::

        sage: from sign_vectors import sign_vector
        sage: sign_vector([5, 0, -1, -2])
        (+0--)
        sage: v = vector([5, 0, -1, -2])
        sage: sign_vector(v)
        (+0--)

    We can also use a string to construct a sign vector::

        sage: sign_vector("00--")
        (00--)
        sage: sign_vector('++-+-00-')
        (++-+-00-)

    Variables are supported to some extent::

        sage: var('a')
        a
        sage: v = vector([1, a, -1])
        sage: sign_vector(v) # not tested TODO fails for some reason
        ...
        UserWarning: Cannot determine sign of symbolic expression, using 0 for sign vector instead.
        (+0-)
        sage: assume(a > 0)
        sage: sign_vector(v)
        (++-)
    """
    if isinstance(iterable, str):
        return SignVector.from_str(iterable)
    support = set()
    psupport = set()
    length = 0
    for entry in iterable:
        sign_entry = sign_symbolic(entry)
        if sign_entry != 0:
            support.add(length)
            if sign_entry > 0:
                psupport.add(length)
        length += 1
    return SignVector(frozenset(support), frozenset(psupport), length)


def zero_sign_vector(length: int):
    r"""
    Return the zero sign vector of a given length.

    INPUT:

    - ``length`` -- length

    EXAMPLES::

        sage: from sign_vectors import zero_sign_vector
        sage: zero_sign_vector(4)
        (0000)
    """
    return SignVector(frozenset(), frozenset(), length)


def random_sign_vector(length: int):
    r"""
    Return a random sign vector of a given length.

    INPUT:

    - ``length`` -- length

    EXAMPLES::

        sage: from sign_vectors import random_sign_vector
        sage: random_sign_vector(5) # random
        (++-0-)

    TEST::

        sage: len(random_sign_vector(5))
        5
    """
    return SignVector.from_str("".join(choices("00+-", k=length)))

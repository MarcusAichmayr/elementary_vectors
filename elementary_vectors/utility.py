"""Utility functions"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################


def is_constant(expression):
    r"""
    Return whether this expression is constant.

    Symbolic expressions are considered constant if they do not depend on any variables.

    EXAMPLES::

        sage: from elementary_vectors.utility import is_constant
        sage: is_constant(5)
        True
        sage: var('a, b')
        (a, b)
        sage: is_constant(a)
        False
        sage: is_constant(-a)
        False
        sage: is_constant(b^2 - a)
        False
        sage: is_constant(SR(5))
        True
    """
    if hasattr(expression, "variables"):
        return not bool(expression.variables())
    return True

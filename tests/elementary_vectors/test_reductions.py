from sage.all import *
from elementary_vectors.reductions import simplify_using_equalities, reduce_factor, reduce_vector, reduce_vectors_support, remove_zero_vectors, reduce_vectors
import unittest

class Tests(unittest.TestCase):
    def setUp(self):
        forget()

    def test_simplify_using_equalities_symbolic(self):
        var('a')
        expr = a + 1
#        self.assertFalse((simplify_using_equalities(expr)).is_constant())

#        assume(a == 0)
        self.assertTrue((simplify_using_equalities(expr, [a == 0])).is_constant())
        
    def test_simplify_using_equalities_polynomial(self):
        R = PolynomialRing(ZZ, 'x')
        x = R.gen()
        assume(SR(x) == 0)
        
        self.assertEqual(simplify_using_equalities(x+1, assumptions()), 1)
        
    def test_reduce_vector(self):
        var('a')
        assume(a == 0)
        v = vector([5*a, 10*a])

        self.assertEqual(reduce_factor(v), vector(v.base_ring(), [1,2]))
        self.assertEqual(reduce_vector(v, eq=assumptions()), zero_vector(v.base_ring(), 2))
        self.assertEqual(reduce_vector(v), vector(v.base_ring(), [1,2]))
        self.assertEqual(reduce_vector(v, eq=[a == 0], factor=False), zero_vector(v.base_ring(), 2))
        self.assertEqual(reduce_vector(v, factor=False), v)
        
    def test_reduce_vectors_support(self):
        l = [vector([1,3,2]), vector([0,0,1]), vector([2,2,0]), vector([0,0,-5])]
        
        self.assertEqual(reduce_vectors_support(l), l[0:3])
        self.assertEqual(reduce_vectors_support([zero_vector(5)]), [zero_vector(5)])

    def test_remove_zero_vectors(self):
        var('a')
        L = [vector([5,2,3]), zero_vector(3), vector([a*0,0,0])]

        self.assertEqual(remove_zero_vectors(L), [vector([5,2,3])])

    def test_reduce_vectors(self):
        var('a')
        assume(a == 0)
        L = [vector([5*a, 10*a, 0]), vector([5*a, 2*a, a]), vector([4, 6, 0])]
        
        self.assertEqual(reduce_vectors(L, eq=assumptions()), [vector([2,3,0])])
        self.assertEqual(reduce_vectors(L), [vector(L[0].base_ring(), [1,2,0]), vector(L[1].base_ring(), [5,2,1])])
        self.assertEqual(reduce_vectors(L, eq=[a == 0], factors=False), [vector([4,6,0])])
        self.assertEqual(reduce_vectors(L, eq=assumptions(), support=True), [vector([2,3,0])])
        self.assertEqual(reduce_vectors(L, eq=assumptions(), support=False, remove_zeros=False), [zero_vector(L[0].base_ring(), 3), zero_vector(L[1].base_ring(),3), vector([2,3,0])])

if __name__ == '__main__':
    unittest.main()

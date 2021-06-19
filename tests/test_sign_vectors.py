from sage.all import *
from sign_vectors import *
import unittest

class SignVectorsTests(unittest.TestCase):
    def setUp(self):
        self.a = sign_vector([1,1,-1,0])
        self.b = sign_vector('+-0+')
        self.c = sign_vector([1,1,0,0])
        self.d = sign_vector('00--')

    def test_getitem_setitem(self):
        e = zero_sign_vector(3)
        self.assertEqual(e[0], 0)
        e[0] = -5 # should set e[0] to -1 = sign(-5)
        self.assertEqual(e[0], -1)
        self.assertEqual(self.a[1], 1) # returns an integer
        self.assertIsInstance(self.a[1:], SignVector) # slicing returns a sign vector

    def test_scalar_multiplication(self):
        self.assertEqual(--self.a, self.a)
        self.assertEqual(1*self.a, self.a)
        self.assertEqual(self.a*(-1), -self.a)
        self.assertEqual(0*self.a, zero_sign_vector(4))

    def test_composition(self):
        self.assertEqual(self.a.compose(self.b), sign_vector([1,1,-1,1]))
        self.assertEqual(self.b.compose(self.a), sign_vector('+--+'))

    def test_conforms(self):
        self.assertEqual(self.a, self.a)
        self.assertNotEqual(self.a, self.b)
        self.assertLess(self.c, self.a)
        self.assertFalse(self.a.conforms(self.b))
        self.assertLess(self.c, self.a)
        self.assertLessEqual(self.c, self.a)
        self.assertGreater(self.a, self.c)
        self.assertGreaterEqual(self.b, self.b)
        
    def test_comparison_with_0(self):
        self.assertEqual(zero_sign_vector(5), 0)
        self.assertLessEqual(zero_sign_vector(5), 0)
        self.assertGreaterEqual(zero_sign_vector(5), 0)
        self.assertLess(self.d, 0)
        self.assertGreater(self.c, 0)
        self.assertNotEqual(self.a, 0)
        self.assertFalse(self.a > 0)
        self.assertFalse(self.a < 0)
        self.assertFalse(0 < self.a)
        # self.assertFalse(self.a.is_zero())
        # self.assertTrue(zero_sign_vector(2).is_zero())
        # self.assertTrue(self.a.is_nonzero())
        with self.assertRaises(TypeError): # same result for real vectors
            self.a > 1
        
    def test_supports(self):
        self.assertEqual(self.a.support(), [0,1,2])
        self.assertEqual(self.a.zero_support(), [3])
        self.assertEqual(self.b.positive_support(), [0,3])
        self.assertEqual(self.b.negative_support(), [1])

    def test_separating_elements(self):
        self.assertEqual(self.a.separating_elements(self.b), [1])
        self.assertEqual(self.c.separating_elements(-self.c), [0,1])

    def test_reverse_signs_in(self):
        self.assertEqual(self.a.reverse_signs_in([0,2]), sign_vector([-1,1,1,0]))
        self.assertEqual(self.a.reverse_signs_in([3]), self.a)

    def test_is_orthogonal(self):
        self.assertTrue(self.a.is_orthogonal_to(self.b))
        self.assertFalse(self.a.is_orthogonal_to(self.c))
        self.assertTrue(self.c.is_orthogonal_to(self.d))
        self.assertTrue(self.a.is_orthogonal_to(zero_sign_vector(4)))

if __name__ == '__main__':
    unittest.main()

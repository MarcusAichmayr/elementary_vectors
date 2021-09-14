from sage.all import *
from sign_vectors import *
from sign_vectors.oriented_matroids import cocircuits_from_matrix
from sign_vectors.utility import loops, is_parallel, parallel_classes, classes_same_support

import unittest

class SignVectorsTests(unittest.TestCase):
    def setUp(self):
        self.a = sign_vector([1,1,-1,0])
        self.b = sign_vector('+-0+')
        T1 = matrix([[0,0,1,-1,0],[1,0,0,0,1],[1,1,1,1,1]])
        T2 = matrix([[0,0,1,-1,0],[1,0,0,0,-1],[1,1,1,1,-1]])
        self.T3 = matrix([[1,-2,0,0,0,0],[2,-4,0,1,0,0],[0,0,1,0,2,3]])
        self.ccT1 = cocircuits_from_matrix(T1)
        self.ccT2 = cocircuits_from_matrix(T2)

    def test_is_parallel_for_sign_vectors(self):
        self.assertTrue(is_parallel(self.ccT1,0,4))
        self.assertFalse(is_parallel(self.ccT1,1,2))
        self.assertFalse(is_parallel(self.ccT1,0,1))
        self.assertTrue(is_parallel(self.ccT2,0,4))

    def test_is_parallel_for_real_vectors(self):
        L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        self.assertTrue(is_parallel(L,0,2))
        self.assertFalse(is_parallel(L,0,1))
        self.assertFalse(is_parallel(L,1,3))
        self.assertTrue(is_parallel(L,4,5))

    def test_is_parallel_ratio(self):
        L = [vector([1,1,2,3,0,0]), vector([-2,1,-4,3,3,-17]), vector([0,1,0,1,0,0])]
        self.assertEqual(is_parallel(L,0,2,return_ratio=True),[True,1/2])
        self.assertEqual(is_parallel(L,0,1,return_ratio=True),[False,0])
        self.assertEqual(is_parallel(L,1,3,return_ratio=True),[False,0])
        self.assertEqual(is_parallel(L,4,5,return_ratio=True),[True,-3/17])

    def test_parallel_classes(self):
        self.assertTrue([0,4] in parallel_classes(self.ccT1))
        self.assertTrue([1] in parallel_classes(self.ccT1))
        self.assertFalse([1,2] in parallel_classes(self.ccT1))
        self.assertTrue([0,1] in parallel_classes(self.T3))
        self.assertTrue([2,4,5] in parallel_classes(self.T3))
        self.assertTrue([3] in parallel_classes(self.T3))

        l = flatten(parallel_classes(self.ccT1))
        l.sort()
        self.assertEqual(l, [0,1,2,3,4]) # each element appears

        self.assertEqual(parallel_classes(self.ccT1), parallel_classes(self.ccT2))

    def test_positive_parallel_classes(self):
        self.assertFalse([0,1] in parallel_classes(self.T3,positive_only=True))
        self.assertTrue([1] in parallel_classes(self.T3,positive_only=True))
        self.assertTrue([0] in parallel_classes(self.T3,positive_only=True))
        
        l = flatten(parallel_classes(self.ccT1,positive_only=True))
        l.sort()
        self.assertEqual(l, [0,1,2,3,4]) # each element appears

    def test_classes_same_support(self):
        self.assertEqual(classes_same_support(self.ccT1)[0][0].support(), classes_same_support(self.ccT1)[0][1].support())
        self.assertNotEqual(classes_same_support(self.ccT1)[2][0].support(), classes_same_support(self.ccT1)[0][1].support())

    def test_loops(self):
        self.assertEqual(loops([sign_vector([0,1,0]), sign_vector([-1,0,0])]), [2])
        self.assertEqual(loops([sign_vector([1,0,0]), sign_vector([-1,0,0])]), [1,2])

    # TODO: do unit tests
    def test_others(self):
        deletion(self.ccT1, [0])
        contraction(self.ccT1, [1])
        closure(self.ccT1)

if __name__ == '__main__':
    unittest.main()

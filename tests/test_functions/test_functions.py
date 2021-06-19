from sage.all import *
from sign_vectors import *
from elementary_vectors import *
from sign_vectors.oriented_matroids import *
import unittest

# tests if appropriate functions are available
class SignVectorsTests(unittest.TestCase):
    def setUp(self):
        self.A = matrix([[1,2,0],[0,1,-1]])
        self.ccA = cocircuits_from_matrix(self.A)

    def test_elementary_vectors(self):
        elementary_vectors(self.A)
        
    def test_oriented_matroids(self):
        covectors_from_cocircuits(self.ccA)
        tA = topes_from_cocircuits(self.ccA)
        cocircuits_from_topes(tA)
        face_enumeration(tA)
        covectors_from_topes(tA)
        topes_from_matrix(self.A)
        topes_from_matrix(self.A, kernel=True)
        covectors_from_matrix(self.A)
        covectors_from_matrix(self.A, kernel=True)
        covectors_from_matrix(self.A, algorithm='fe')
        covectors_from_matrix(self.A, algorithm='fe', separate=True)

    def test_others(self):
        deletion(self.ccA, [0])
        contraction(self.ccA, [1])
        conformal_elimination(vector([1,2]), vector([-1,3]))
        closure(self.ccA)

if __name__ == '__main__':
    unittest.main()

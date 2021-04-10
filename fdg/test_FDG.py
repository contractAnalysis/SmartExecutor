from fdg.FDG import FDG
from fdg.sequence import Sequence
import unittest

class TestSequence(unittest.TestCase):

    def test_get_combination_sequences(self):
        data=[[1,4], [2,6],[5]]
        seq=Sequence()
        result=seq.get_combination(data)
        expected=[(1, 2), (1, 6), (4, 2), (4, 6), (1, 5), (4, 5), (2, 5), (6, 5), (1, 2, 5), (1, 6, 5), (4, 2, 5), (4, 6, 5)]
        self.assertEqual(result,expected,"does not equal")

if __name__=='__main__':

    unittest.main()

from fdg.FDG import FDG
from fdg.sequence import Sequence
import unittest

class TestSequence(unittest.TestCase):

    def test_get_parents_combination(self):
        functions_dict = {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'],
                          'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'],
                          'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'],
                                 '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'],
                          'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'],
                          'f6': ['decreaseApproval', [], ['allowed'], '66188463'],
                          'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'],
                          'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'],
                                 '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'],
                          'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'],
                          'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'],
                          'f0': ['constructor', [], ['owner'], '8afc3605']}

        fdg = FDG(functions_dict)
        seq = Sequence(fdg, 'burn')
        result=seq.get_parents_combination()
        expected={2: [(3, 10), (8, 10)]}
        self.assertEqual(result, expected, "does not equal")

    def test_dict_ftns_short_sequences(self):
        functions_dict = {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'],
                          'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'],
                          'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'],
                                 '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'],
                          'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'],
                          'f6': ['decreaseApproval', [], ['allowed'], '66188463'],
                          'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'],
                          'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'],
                                 '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'],
                          'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'],
                          'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'],
                          'f0': ['constructor', [], ['owner'], '8afc3605']}

        fdg = FDG(functions_dict)
        seq = Sequence(fdg, 'burn')
        re=seq.dict_ftns_short_sequences()
        expected = {1: {10: [[10]]}, 2: {3: [[5, 3], [6, 3]], 8: [[7, 8]], 10: [[1, 10]]}, 3: {2: [[5, 3, 2], [6, 3, 2], [1, 7, 8, 2], [1, 7, 9, 2]], 3: [[5, 3], [6, 3]], 8: [[1, 7, 8]]}}


        self.assertEqual(re, expected, "does not equal")

    def test_generate_depth_seqeunces(self):
        functions_dict = {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'],
                          'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'],
                          'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'],
                                 '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'],
                          'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'],
                          'f6': ['decreaseApproval', [], ['allowed'], '66188463'],
                          'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'],
                          'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'],
                                 '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'],
                          'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'],
                          'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'],
                          'f0': ['constructor', [], ['owner'], '8afc3605']}

        fdg = FDG(functions_dict)
        seq = Sequence(fdg, 'burn')
        result=seq.generate_depth_seqeunces()
        expected={2: [[3, 10, 11], [8, 10, 11]]}
        self.assertEqual(result, expected, "does not equal")

if __name__=='__main__':

    unittest.main()





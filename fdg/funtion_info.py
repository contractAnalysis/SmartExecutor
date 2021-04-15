
from slither.slither import Slither
from slither.utils.function import get_function_id
import sha3
class Function_info():
    def __init__(self,solidity_file:str, contract_name:str):
        self.solidity_file=solidity_file
        self.contract_name=contract_name

    def functions_dict_slither(self):
        '''
        get function info
        :param contract_file:
        :param contract_name:
        :return: list of (function: r_sv_list: w_sv_list)
        '''
        # Init slither
        slither = Slither(self.solidity_file)

        # Get the contract
        contract = slither.get_contract_from_name(self.contract_name)
        assert contract

        constructors_list = []
        for item in contract.constructors:
            constructors_list.append(item.name)

        function_dict = {}

        i = 1
        for f in contract.functions:
            if f.name.__eq__('slitherConstructorVariables'): continue
            if f.name.__eq__('slitherConstructorConstantVariables'): continue
            # if f.name in constructors_list:
            #     print(f'{f.name} belongs to part of constructor!')
            #     continue
            if f.is_constructor: continue



            func_hash = self.get_function_id(f.full_name)

            r_list = []
            f_r = f.all_conditional_state_variables_read()
            if len(f_r) > 0:
                r_list = [sv.name for sv in f_r]

            w_list = []
            f_w = f.all_state_variables_written()
            if len(f_w) > 0:
                w_list = [sv.name for sv in f_w]

            if len(r_list) != 0 or len(w_list) != 0:
                function_dict['f' + str(i)] = [f.full_name, r_list, w_list, func_hash]
                i = i + 1

        # check if constructor is included or not. if not, add it
        constructor_write_list = []
        f_w = []
        for item in contract.constructors:
            f_w += item.all_state_variables_written()
        if len(f_w) > 0:
            constructor_write_list = [sv.name for sv in f_w]

        if len(contract.constructors)>0:
            func_hash = str(hex(get_function_id(contract.constructor.full_name))).strip('0x')
            function_dict['f0'] = ["constructor", [], constructor_write_list,
                               func_hash]  # abstract constructor: a combination of specific constructors from different contracts with inheritant relationship

        else: # in the case of no constructor
            function_dict['f0']=["constructor",[],[],""]


        return function_dict

    def get_function_id(self,sig: str) ->str:
        """'
            Return the function id of the given signature
        Args:
            sig (str)
        Return:
            (int)
        """
        s = sha3.keccak_256()
        s.update(sig.encode("utf-8"))
        return s.hexdigest()[:8]

if __name__=='__main__':
    # ftn_info=Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/wei_test.sol', 'wei_test')

    ftn_dict= ftn_info.functions_dict_slither()
    print("===== ftn_dict ====")
    for key, value in ftn_dict.items():
        print("\t{}:  {}".format(key, value))
    pass


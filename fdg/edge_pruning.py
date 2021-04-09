

from slither.core.expressions import UnaryOperation, BinaryOperation, CallExpression, Identifier, Literal, \
    AssignmentOperation, TupleExpression, BinaryOperationType
from slither.core.solidity_types import ElementaryType
from slither.slither import Slither

from slither.utils.function import get_function_id
import z3

from z3 import Solver

import numpy as np
def try_(contract_file:str, contract_name:str):
    slither = Slither(contract_file)
    # Get the contract
    contract = slither.get_contract_from_name(contract_name)
    assert contract

    sv_initialized={}
    # mapping,array,have no initial values
    sv_constant=[]
    for sv in contract.state_variables:
        if sv.is_constant:
            sv_constant.append(sv.expression.__str__())
        if sv.initialized:
            if isinstance(sv.expression,BinaryOperation):
                # calculate the final result and save
                nested_exp_list = visit_exp_recurcively(sv.expression)
                if sv.type.__str__().__contains__('int'):
                    result = calculate_from_nested_list(nested_exp_list,'int')
                else:
                    result = calculate_from_nested_list(nested_exp_list,'float')
                sv_initialized[sv.name.__str__()] = result
            elif isinstance(sv.expression,Literal):
                sv_initialized[sv.name.__str__()]=sv.expression.__str__()
            else: pass
        # assume default value for integer variable
        elif isinstance(sv.type,ElementaryType):
            if sv.type.__str__().__contains__('int'):
                # assume that the default value for integer variable is 0
                sv_initialized[sv.name]='0'


    # consider the initialization in constructor
    for item in contract.constructors:
        for exp in item.expressions:
            if isinstance(exp,AssignmentOperation):
                # now only consider integer variable
                if exp.expressions[0].value.type.__str__().__contains__('int'):
                    nested_exp=visit_exp_recurcively(exp.expressions[1])
                    re=calculate_from_nested_list_1(nested_exp,'int')
                    sv_name=exp.expressions[0].value
                    sv_initialized[exp.expressions[0].value.__str__()]=re
                # print(f'{re}')
                # print(f'{exp.expressions[0].value}')

    print('==== state variables that have initial values ====')
    for key,value in sv_initialized.items():
        print(f'\t{key}: {value}')
    print()

    print('==== constant state variables ====')
    for item in sv_constant:
        print(f'\t {item}')
    print()

    for ftn in contract.functions:
        print(f'==== function: {ftn.name} ====')
        r_list = []
        f_r = ftn.all_conditional_state_variables_read()
        if len(f_r) > 0:
            r_list = [sv.name for sv in f_r]
        for item in f_r:
            print(f'\t sv read:  {item}')


        # get expressions containing state variables that have initial values
        target_exp=[]
        # get expressions from requires
        for exp in ftn.expressions:
            if isinstance(exp,CallExpression):
                for sv in f_r:
                    if sv.__str__() in sv_initialized.keys():
                        if exp.__str__().__contains__(sv.__str__()):
                            if exp not in target_exp:
                                target_exp.append(exp)
            # capture conditions in if, for, while statements ( currently, To-Do)
            elif isinstance(exp,BinaryOperation):
                for sv in f_r:
                    if exp.__str__().__contains__(sv.__str__()):
                        target_exp.append(exp)

            else: pass
        # get expressions from modifiers
        for modifier in ftn.modifiers:
            for exp in modifier.expressions:
                if isinstance(exp,CallExpression):
                    for sv in f_r:
                        if sv.__str__() in sv_initialized.keys():
                            if exp.__str__().__contains__(sv.__str__()):
                                if exp not in target_exp:
                                    target_exp.append(exp)
                # capture conditions in if, for, while statements ( currently, To-Do)
                elif isinstance(exp,BinaryOperation):
                    for sv in f_r:
                        if exp.__str__().__contains__(sv.__str__()):
                            target_exp.append(exp)

                else: pass

        conditions_list=[]
        if len(target_exp)!=0:
            for t_exp in target_exp:
                print(f'\ttarget_exp:  {t_exp}')
                for arg in t_exp.arguments:
                    re= visit_exp_recurcively(arg)

                    cons=[]
                    get_conditions_from_nested_list(re,cons)
                    # print(cons)
                    conditions_list+=cons


        if len(conditions_list)!=0:
            print(f'\tbasic conditions:')
            for item in conditions_list:
                if len(item)!=0:
                    print(f'\t\t{item}')


        sv_=get_sv_sat_conditions(conditions_list,sv_initialized,sv_constant)
        if len(sv_)!=0:
            print(f'\tstate variables that can be removed:')
            print(f'\t\t{sv_}')

        final_f_r=list(set(r_list)-set(sv_))

        if len(r_list)!=0:
            print(f'\t++++++++++++++++++++++++++++++')
            print(f'\t depends on  previously:')
            for item in r_list:
                print(f'\t\t{item}')

            print(f'\t depends on now:')
            if len(final_f_r)==0:
                print(f'\t\t{[]}')
            else:
                for item in final_f_r:
                    print(f'\t\t{item}')


#recurcively to get targets in a nested list. a target is a condition.
# get a list of individual conditions
def get_conditions_from_nested_list(nested_list,result:list):
    if isinstance(nested_list, list):
        # in the case of BinaryOperation: len==3,depth==2
        if len(nested_list)==3 and depth_of_list()(nested_list) == 2:
            result.append(nested_list)

        # in the case of UnaryOperation: len==2,depth==2 ( a boolean variable and a unary operation)
        elif len(nested_list)==2 and depth_of_list()(nested_list)==2:
            result.append(nested_list)

        # in the case of UnaryOperation:len==1,depth==1 (only a boolean variable as condition)
        elif len(nested_list)==1 and depth_of_list()(nested_list)==1:
            result.append(nested_list)
        else:
            if len(nested_list)==1 and depth_of_list()(nested_list)>1:
                nested_list=nested_list[0]
                get_conditions_from_nested_list(nested_list,result)
            else:
                for item in nested_list:
                    get_conditions_from_nested_list(item,result )

            # for item in nested_list:
            #     show_re(item)

def get_sv_sat_conditions(condition_list,sv_initialized:dict,sv_constant:list):
    sv_list=[]
    # help to remove sv that satisfies one condition but does not in another condition
    sv_unsat_list=[]
    for cond in condition_list:
        if len(cond)!=0:
            # format: ['type: EQUAL', ['None: phase'], ['uint256: 2']]
            status,sv=sat_check(cond,sv_initialized)
            if str(status)=='sat':
                sv_list+=sv
            else:
                sv_unsat_list+=sv
    return list((set(sv_list)-set(sv_unsat_list)).union(set(sv_constant)))

#recurcively to print targets in a nested list
def show_re(re):
    if isinstance(re, list):
        # in the case of BinaryOperation: len==3,depth==2
        if len(re)==3 and depth_of_list()(re) == 2:
            print(f'\t{re}')
        # in the case of UnaryOperation: len==2,depth==2 ( a boolean variable and a unary operation)
        elif len(re)==2 and depth_of_list()(re)==2:
            print(f'\t{re}')
        # in the case of UnaryOperation:len==1,depth==1 (only a boolean variable as condition)
        elif len(re)==1 and depth_of_list()(re)==1:
            print(f'\t{re}')
        else:
            for item in re:
                show_re(item)
    # else:
    #     if isinstance(re, str):
    #         print('this is a string')


def depth_of_list():
    # get the depth of a nested list
    a = lambda x: isinstance(x, list) and max(map(a, x)) + 1
    return a

#recurcively to visit an expressions and return the results as a nested list
def visit_exp_recurcively(exp):
    re=[]
    if isinstance(exp,Identifier):
        # re.append(f'{exp.value.type}: {exp.value}')
        re.append(f'{exp.type}:{exp.value}:{exp.value.type}')
    elif isinstance(exp,Literal):
        re.append(f'{exp.type}:{exp.value}')
    else:
        if isinstance(exp,BinaryOperation):
            re.append(f'type:{exp.type.name}')
            for exp1 in exp.expressions:
                re.append(visit_exp_recurcively(exp1))

        elif isinstance(exp,UnaryOperation):
            re.append(f'type:{exp.type.name}')
            re.append(visit_exp_recurcively(exp.expression))

        elif isinstance(exp,TupleExpression):
            for exp3 in exp.expressions:
                re.append(visit_exp_recurcively(exp3))
        else:
            re.append(print(f'type:{exp.type.name}'))
            for exp4 in exp.expressions:
                re.append(visit_exp_recurcively(exp4))

    return re

#recurcively to print expressions
def print_exp_0(exp,i):

    print(f'i={i}')
    i+=1
    if isinstance(exp,Identifier) or isinstance(exp,Literal):
        print(f'\t Identifier or Literal: {exp.value}')
    else:

        if isinstance(exp,BinaryOperation):
            print(f'\t type: {exp.type.name}')
            for exp1 in exp.expressions:
                print_exp_0(exp1,i)
        elif isinstance(exp,UnaryOperation):
            print(f'\t type:{exp.type.name}')
            if isinstance(exp.expression, Identifier) or isinstance(exp.expression, Literal):
                print(f'\t Identifier or Literal: {exp.expression.value}')
            else:
                for exp2 in exp.expression:
                    print_exp_0(exp2,i)

        elif isinstance(exp,TupleExpression):
            for exp3 in exp.expressions:
                print_exp_0(exp3,i)
        else:
            print(f'\t type: {exp.type.name}')
            for exp4 in exp.expressions:
                print_exp_0(exp4,i)








def get_sv_should_be_removed(contract_file:str, contract_name:str):
    slither = Slither(contract_file)
    # Get the contract
    contract = slither.get_contract_from_name(contract_name)
    assert contract

    sv_initialized={}
    # mapping,array,have no initial values
    sv_constant=[]
    for sv in contract.state_variables:
        if sv.is_constant:
            sv_constant.append(sv.expression.__str__())
        if sv.initialized:
            if isinstance(sv.expression,BinaryOperation):
                # calculate the final result and save
                nested_exp_list = visit_exp_recurcively(sv.expression)
                if sv.type.__str__().__contains__('int'):
                    result = calculate_from_nested_list(nested_exp_list,'int')
                else:
                    result = calculate_from_nested_list(nested_exp_list,'float')
                sv_initialized[sv.name.__str__()] = result
            elif isinstance(sv.expression,Literal):
                sv_initialized[sv.name.__str__()]=sv.expression.__str__()
            else: pass
        # assume default value for integer variable
        elif isinstance(sv.type,ElementaryType):
            if sv.type.__str__().__contains__('int'):
                # assume that the default value for integer variable is 0
                sv_initialized[sv.name]='0'


    # consider the initialization in constructor
    for item in contract.constructors:
        for exp in item.expressions:
            if isinstance(exp,AssignmentOperation):
                # now only consider integer variable
                if exp.expressions[0].value.type.__str__().__contains__('int'):
                    nested_exp=visit_exp_recurcively(exp.expressions[1])
                    re=calculate_from_nested_list_1(nested_exp,'int')
                    sv_name=exp.expressions[0].value
                    sv_initialized[exp.expressions[0].value.__str__()]=re
                # print(f'{re}')
                # print(f'{exp.expressions[0].value}')

    ftn_sv_read_removed={} #save state variables that shuould be removed for each function
    for ftn in contract.functions:
        r_list = []
        f_r = ftn.all_conditional_state_variables_read()
        if len(f_r) > 0:
            r_list = [sv.name for sv in f_r]

        # get expressions containing state variables that have initial values
        target_exp=[]
        # get expressions from requires
        for exp in ftn.expressions:
            if isinstance(exp,CallExpression):
                for sv in f_r:
                    if sv.__str__() in sv_initialized.keys():
                        if exp.__str__().__contains__(sv.__str__()):
                            if exp not in target_exp:
                                target_exp.append(exp)
            # capture conditions in if, for, while statements ( currently, To-Do)
            elif isinstance(exp,BinaryOperation):
                for sv in f_r:
                    if exp.__str__().__contains__(sv.__str__()):
                        target_exp.append(exp)

            else: pass
        # get expressions from modifiers
        for modifier in ftn.modifiers:
            for exp in modifier.expressions:
                if isinstance(exp,CallExpression):
                    for sv in f_r:
                        if sv.__str__() in sv_initialized.keys():
                            if exp.__str__().__contains__(sv.__str__()):
                                if exp not in target_exp:
                                    target_exp.append(exp)
                # capture conditions in if, for, while statements ( currently, To-Do)
                elif isinstance(exp,BinaryOperation):
                    for sv in f_r:
                        if exp.__str__().__contains__(sv.__str__()):
                            target_exp.append(exp)
                else: pass

        # recursively visit expressions and get conditions
        conditions_list=[]
        if len(target_exp)!=0:
            for t_exp in target_exp:
                for arg in t_exp.arguments:
                    re= visit_exp_recurcively(arg)
                    cons=[]
                    get_conditions_from_nested_list(re,cons)
                    conditions_list+=cons
        # z3 checks satisfiability to get state variables that makes conditions satisfied
        sv_=get_sv_sat_conditions(conditions_list,sv_initialized,sv_constant)

        if len(sv_)!=0:
            ftn_sv_read_removed[ftn.name]=sv_
    return ftn_sv_read_removed

def functions_dict_slither_prune(functions_dict:dict,functions_sv_should_be_removed:dict):
    for key,value in functions_dict.items():
        if isinstance(value,list):
            if value[0] in functions_sv_should_be_removed.keys():
                functions_dict.get(key)[1]=list(set(functions_dict.get(key)[1])-set(functions_sv_should_be_removed.get(value[0])))



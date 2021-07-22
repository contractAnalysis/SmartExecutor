import sha3
def get_function_id(sig: str) -> str:
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


# get valid pc interval from a list of PCs of opcode GT
def get_valid_pc_interval(gt_pc_list:list,max_pc_value:int):
    gt_pc_list.sort()
    pairs = []
    step = 5

    for i in range(1, len(gt_pc_list)):
        if gt_pc_list[i] == gt_pc_list[i - 1] + step:
            continue
        else:
            pairs.append((gt_pc_list[i - 1], gt_pc_list[i]))
    pairs.append((gt_pc_list[-1], max_pc_value))
    return pairs

# check if a pc is in the valid pc intervals
def pc_is_valid(pc:int,valid_pc_intervals:list):
    for item in valid_pc_intervals:
        if pc >item[0] and pc<item[1]:
            return True
    return False

def assign_pc_seq_exe_phase_dup1(current_pc: int, pc_list: list, pc_interval_dict: dict) -> int:
    if len(pc_list)>0:
        for pc in pc_list:
            if current_pc <= pc:
                return pc
    return pc_interval_dict['pc_interval_end']

def assign_pc_fdg_phase_dup1(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict, depth:int,ftn_not_covered:dict,pc_interval_dict: dict) -> int:
    if prt_ftn_name in ftn_to_index_dict.keys():
        prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
        if prt_ftn_idx in fdg_pc_dict.keys():
            pc_list=fdg_pc_dict[prt_ftn_idx]
            if len(pc_list)>0:
                for pc in pc_list:
                    if current_pc <= pc:
                        return pc
                # return pc_interval_dict['pc_interval_end']

    # fallback, no-edge nodes, leaf nodes
    if depth in ftn_not_covered.keys():
        ftn_uncovered=ftn_not_covered[depth]
        if len(ftn_uncovered)>=1:
            for pc in ftn_uncovered:
                if current_pc<=pc:
                    return pc

    if 'fallback' in ftn_not_covered.keys():
        for pc in ftn_not_covered['fallback']:
            if current_pc<=pc:
                return pc
    return pc_interval_dict['pc_interval_end']






# def assign_pc_fdg_phase_jumpi(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict, depth:int,ftn_not_covered:dict,pc_interval_dict: dict) -> int:
#     if prt_ftn_name in ftn_to_index_dict.keys():
#         prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
#         if prt_ftn_idx in fdg_pc_dict.keys():
#             pc_list=fdg_pc_dict[prt_ftn_idx]
#             if len(pc_list)>=1:
#                 temp3=pc_list+[pc_interval_dict['pc_interval_end']]
#                 for pc in temp3:
#                     if current_pc < pc:
#                         return pc
#                 return current_pc
#
#     # when prt_ftn_name does not in FDG
#     # lead nodes
#     if depth in ftn_not_covered.keys():
#         ftn_uncovered=ftn_not_covered[depth]
#         if len(ftn_uncovered)>=1:
#             temp=ftn_uncovered+[pc_interval_dict['pc_interval_end']]
#             for pc in temp:
#                 if current_pc<pc:
#                     return pc
#             return current_pc
#     return pc_interval_dict['pc_interval_end']
#
# def assign_pc_fdg_phase_dup1(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict, depth:int,ftn_not_covered:dict,pc_interval_dict: dict) -> int:
#     if prt_ftn_name in ftn_to_index_dict.keys():
#         prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
#         if prt_ftn_idx in fdg_pc_dict.keys():
#             pc_list=fdg_pc_dict[prt_ftn_idx]
#             if len(pc_list)>=1:
#                 return pc_list[0]
#
#
#     # when prt_ftn_name does not in FDG
#     # handle leaf nodes
#     if depth in ftn_not_covered.keys():
#         ftn_uncovered=ftn_not_covered[depth]
#         if len(ftn_uncovered)>0:
#             return ftn_uncovered[0]
#         else:return pc_interval_dict['pc_interval_end']
#
#     return pc_interval_dict['pc_interval_end']
#



# def assign_pc_seq_exe_phase_jumpi(current_pc: int, pc_list: list, pc_interval_dict: dict,ftn_wo_edges_not_covered_pc:list) -> int:
#
#     temp1=pc_list+ftn_wo_edges_not_covered_pc
#     temp1.sort()
#     if len(temp1)==0:
#         return pc_interval_dict['pc_interval_end']
#     else:
#         temp1+=[pc_interval_dict['pc_interval_end']]
#         for pc in temp1:
#             if current_pc < pc:
#                 return pc
#         return current_pc
#
# def assign_pc_seq_exe_phase_dup1(current_pc: int, pc_list: list, pc_interval_dict: dict,ftn_wo_edges_not_covered_pc:list) -> int:
#
#     temp2 = pc_list + ftn_wo_edges_not_covered_pc
#     temp2.sort()
#
#     if len(temp2)==0:
#         return pc_interval_dict['pc_interval_end']
#     else:
#         return temp2[0]
#
# def assign_pc_fdg_phase_jumpi(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict, depth:int,ftn_not_covered:dict,pc_interval_dict: dict,ftn_wo_edges_not_covered_pc:list) -> int:
#     if prt_ftn_name in ftn_to_index_dict.keys():
#         prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
#         if prt_ftn_idx in fdg_pc_dict.keys():
#             pc_list=fdg_pc_dict[prt_ftn_idx]
#             temp3 = pc_list + ftn_wo_edges_not_covered_pc
#             temp3.sort()
#
#             if len(temp3)>=1:
#                 temp3+=[pc_interval_dict['pc_interval_end']]
#                 for pc in temp3:
#                     if current_pc < pc:
#                         return pc
#                 return current_pc
#         else: # no dependent children
#             if len(ftn_wo_edges_not_covered_pc)>0:
#                 temp_=ftn_wo_edges_not_covered_pc+[pc_interval_dict['pc_interval_end']]
#                 for pc in temp_:
#                     if current_pc < pc:
#                         return pc
#                 return current_pc
#             else: return pc_interval_dict['pc_interval_end']
#
#         return pc_interval_dict['pc_interval_end']
#     # when prt_ftn_name does not in FDG
#     # or it does not have edges(child nodes)
#     if depth in ftn_not_covered.keys():
#         ftn_uncovered=ftn_not_covered[depth]
#         if len(ftn_uncovered)>=1:
#             temp=ftn_uncovered+[pc_interval_dict['pc_interval_end']]
#             for pc in temp:
#                 if current_pc<pc:
#                     return pc
#             return current_pc
#
#     return pc_interval_dict['pc_interval_end']
#
# def assign_pc_fdg_phase_dup1(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict, depth:int,ftn_not_covered:dict,pc_interval_dict: dict,ftn_wo_edges_not_covered_pc:list) -> int:
#     if prt_ftn_name in ftn_to_index_dict.keys():
#         prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
#         if prt_ftn_idx in fdg_pc_dict.keys():
#             pc_list=fdg_pc_dict[prt_ftn_idx]
#
#             temp4 = pc_list + ftn_wo_edges_not_covered_pc
#             temp4.sort()
#
#             if len(temp4)>=1:
#                 return temp4[0]
#             else: return pc_interval_dict['pc_interval_end']
#         else: # no dependent children
#             if len(ftn_wo_edges_not_covered_pc)>0:
#                 return ftn_wo_edges_not_covered_pc[0]
#             else: return pc_interval_dict['pc_interval_end']
#
#
#     # when prt_ftn_name does not in FDG
#     if depth in ftn_not_covered.keys():
#         ftn_uncovered=ftn_not_covered[depth]
#         if len(ftn_uncovered)>0:
#             return ftn_uncovered[0]
#         else:return pc_interval_dict['pc_interval_end']
#
#     return pc_interval_dict['pc_interval_end']
#


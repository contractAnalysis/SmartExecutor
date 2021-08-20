
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

def assign_pc_seq_exe_phase_dup1(current_pc: int, pc_list: list,ftn_not_covered_pc:list, pc_interval_dict: dict) -> int:
    if len(pc_list)>0:
        for pc in pc_list:
            if current_pc <= pc:
                return pc
        return pc_interval_dict['pc_interval_end']
    if len(ftn_not_covered_pc)>=1:
        for pc in ftn_not_covered_pc:
            if current_pc<=pc:
                return pc
    return pc_interval_dict['pc_interval_end']

def assign_pc_fdg_phase_dup1(current_pc: int, prt_ftn_name:str,ftn_to_index_dict:dict,fdg_pc_dict:dict,ftn_not_covered_pc:list,pc_interval_dict: dict) -> int:
    if prt_ftn_name in ftn_to_index_dict.keys():
        prt_ftn_idx=ftn_to_index_dict[prt_ftn_name]
        if prt_ftn_idx in fdg_pc_dict.keys():
            pc_list=fdg_pc_dict[prt_ftn_idx]
            if len(pc_list)>0:
                for pc in pc_list:
                    if current_pc <= pc:
                        return pc
                return pc_interval_dict['pc_interval_end']

    # when prt_ftn_name does not in FDG
    # lead nodes

    if len(ftn_not_covered_pc)>=1:
        for pc in ftn_not_covered_pc:
            if current_pc<=pc:
                return pc
    return pc_interval_dict['pc_interval_end']




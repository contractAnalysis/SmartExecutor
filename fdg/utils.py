import itertools as it

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


def assign_pc(current_pc: int, pc_list: list,executed_pc:list,ftn_not_covered_pc:list, pc_interval_dict: dict) -> int:
    if len(pc_list)>0:
        for pc in pc_list:
            if pc in executed_pc:continue
            if current_pc <= pc:
                return pc
        return pc_interval_dict['pc_interval_end']
    if len(ftn_not_covered_pc)>=1:
        for pc in ftn_not_covered_pc:
            if pc in executed_pc: continue
            if current_pc<=pc:
                return pc
    return pc_interval_dict['pc_interval_end']

def assign_pc_fdg_phase_dup1(current_pc: int, pc_list:list,pc_interval_dict: dict) -> int:
    if len(pc_list)>0:
        for pc in pc_list:
            if current_pc <= pc:
                return pc
    return pc_interval_dict['pc_interval_end']

def assign_pc_special_ftn_dup1(current_pc: int, pc_list:list,executed_pc:list,pc_interval_dict: dict) -> int:
    if len(pc_list)>0:
        for pc in pc_list:
            if pc in executed_pc:continue
            if current_pc <= pc:
                return pc
    return pc_interval_dict['pc_interval_end']



def get_combination(list_for_comb,comb_length:int):
    """

    :param list_for_comb: [[1,4], [2,6],[5]]
    :param comb_length: 2 (two elements in a combination)
    :return: [(1, 2), (1, 6), (4, 2), (4, 6), (1, 5), (4, 5), (2, 5), (6, 5)]
    """
    com_re = []
    # do combination with length
    num_groups = len(list_for_comb)
    if num_groups<comb_length:return []

    # get group combinations
    com_groups = it.combinations(list_for_comb, comb_length)

    for groups in com_groups:
        com_re +=it.product(*list(groups))

    return com_re


# merge two lists, common elements are considered once
# the order in each list is maintained in combined list
# the first element in the combined result should be the first element in a



def merge_fix_list_1_specified_lenth(list_1, list_2, specified_length):
    if specified_length >= len(list_1):
        specified_length = len(list_1)
    c = list_1[0:specified_length]  # get the fixed prefix
    j = 0
    # remove elements of list_2 appearing in fixed prefix
    while j < len(list_2):
        for i in range(specified_length):
            if j == len(list_2):
                break
            if list_1[i] == list_2[j]:
                j += 1
        break
    a = list_1[specified_length:]  # get elements after the fixed prefix
    b = list_2[j:]  # get the left behind elements in list_2
    re = merge_two_list(a, b)  # merge two list
    if re:
        c += re
    return c

def merge_fix_list_1_specified_lenth_no_repeated_elements(list_1, list_2, specified_length):
    if specified_length >= len(list_1):
        specified_length = len(list_1)
    c = list_1[0:specified_length]  # get the fixed prefix
    j = 0
    # remove elements of list_2 appearing in fixed prefix
    list_2_left=[]
    while j<len(list_2):
        if list_2[j] not in list_1[0:specified_length]:
           list_2_left.append(list_2[j])
        j+=1
    list_1_left = list_1[specified_length:]  # get elements after the fixed prefix
    re = merge_two_list(list_1_left, list_2_left)  # merge two list
    if re:
        c += re
    return c

def merge_two_list(list_1, list_2) -> list:
    c = []  # save the merged result
    flag_common_element=False
    while len(list_1)>0 and len(list_2)>0:
        for i in range(len(list_1)):
            for j in range(len(list_2)):
                if list_1[i]==list_2[j]:
                    flag_common_element=True
                    break
            if flag_common_element:
                break
        if flag_common_element:
            flag_common_element=False  # add this statement to remove the bug in this method
            c+=list_1[0:i]
            c+=list_2[0:j+1]
            list_1=list_1[i+1:]
            list_2=list_2[j+1:]
        else:
            c += list_1
            c += list_2
            return c
    if len(list_1)>0:
        c+=list_1
    if len(list_2)>0:
        c+=list_2
    return c



if __name__ == '__main__':
    # print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [5], 3)}')
    #
    # # #print(merge_fix_list_1([1,2,5,7],[4,5]))
    # # #print(merge_fix_list_1([1, 2, 5, 7], [2,8]))
    # # #print(merge_fix_list_1([1, 2, 5, 7], [2,7]))
    # # #print(merge_fix_list_1([1, 2, 5, 7], [2,1]))
    # # #print(merge_fix_list_1([1, 2, 5, 7], [1,3,7]))
    # # print('=====================')
    #
    # # print(f'{[1,2,5,7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1, 5],3)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],5)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],4)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],3)}')
    # print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 2)}')
    # print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 1)}')
    # print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 0)}')
    #
    # print('=====================')
    # print(f'{[1, 2, 5, 7, 3, 5, 6, 9]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5, 6, 9], 5)}')
    # print(f'{[1, 2, 5, 7, 3, 5, 6, 9]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5, 6, 9], 4)}')
    # print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [5], 3)}')
    # print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1], 2)}')
    # print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [], 2)}')
    # print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 0)}')

    pass





    print(f'{[7,3,5]}--{merge_two_list([7],[3,5])}')
    print(f'{[3,5,7]}--{merge_two_list([3,5],[7])}')
    print(f'{[3,5]}--{merge_two_list([3,5],[3])}')
    print(f'{[3,5]}--{merge_two_list([3,5],[5])}')
    print(f'{[3,5]}--{merge_two_list([3],[3,5])}')
    print(f'{[3,5]}--{merge_two_list([5],[3,5])}')
    print(f'{[3,5]}--{merge_two_list([3,5],[3,5])}')
    print(f'{[3,5]}--{merge_two_list([],[3,5])}')
    print(f'{[3,5]}--{merge_two_list([3,5],[])}')
    print(f'{[]}--{merge_two_list([],[])}')

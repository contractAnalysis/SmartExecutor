
# merge two lists, common elements are considered once
# the order in each list is maintained in combined list
# the first element in the combined result should be the first element in a
def merge_2_list(a, b) -> list:
    a1 = []  # temporarily save elements in a that are not matched with element in b
    b1 = []  # temporarily save elements in b that are not matched with element in a
    c = []  # save the combined result

    if len(b) == 0:
        c = a
        return c
    if len(a) == 0:
        c = b
        return c
    j = 0
    while j < len(b):
        flag = True
        for i in range(len(a)):
            if b[j] == a[i]:
                # add elements to c
                if len(a1) > 0:
                    c += a1
                    a1 = []
                if len(b1) > 0:
                    c += b1
                    b1 = []
                c.append(a[i])
                # update a
                a = a[i + 1:]  # remove elements that already moved to c
                # check if either a or b reaches end
                if len(a) == 0:
                    if j < len(b) - 1:
                        c += b[j + 1:]
                    return c
                else:
                    if j == len(b) - 1:
                        c += a
                        return c
                    else:  # both a and b still do not reach end, go to the next iteration
                        j += 1
                flag = False
                break
            else:  # do not match, save a[i] to a1
                a1.append(a[i])

        if flag:
            if j < len(b) - 1:
                a1 = []
                b1.append(b[j])
                j += 1
            else:
                c += a
                c += b1
                c.append(b[j])
                return c


def merge_fix_list_1(list_1, list_2):
    c = list_1
    j = 0
    # remove elements of list_2 that appear in list_1
    while j < len(list_2):
        for i in range(len(list_1)):
            if j == len(list_2): break
            if list_1[i] == list_2[j]:
                j += 1

        if j < len(list_2):
            c += list_2[j:]
        break
    return c


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
    re = merge_2_list(a, b)  # merge two list
    if re:
        c += re
    return c


if __name__ == '__main__':
    print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [5], 3)}')

    # #print(merge_fix_list_1([1,2,5,7],[4,5]))
    # #print(merge_fix_list_1([1, 2, 5, 7], [2,8]))
    # #print(merge_fix_list_1([1, 2, 5, 7], [2,7]))
    # #print(merge_fix_list_1([1, 2, 5, 7], [2,1]))
    # #print(merge_fix_list_1([1, 2, 5, 7], [1,3,7]))
    # print('=====================')

    # # print(f'{[1,2,5,7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1, 5],3)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],5)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],4)}')
    # print(f'{[1,2,5,7,3,5]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7],[1,3,5],3)}')
    print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 2)}')
    print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 1)}')
    print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 0)}')

    print('=====================')
    print(f'{[1, 2, 5, 7, 3, 5, 6, 9]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5, 6, 9], 5)}')
    print(f'{[1, 2, 5, 7, 3, 5, 6, 9]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5, 6, 9], 4)}')
    print(f'{[1, 2, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [5], 3)}')
    print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1], 2)}')
    print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [], 2)}')
    print(f'{[1, 2, 3, 5, 7]}--{merge_fix_list_1_specified_lenth([1, 2, 5, 7], [1, 3, 5], 0)}')

    print('=====================')

    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4, 5, 6, 9], [4, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4, 5, 6, 9], [4, 6])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4, 5, 6, 9], [7, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4, 5, 6, 9], [9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4, 5, 6, 9], [7])}')
    print('=====================')

    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 4], [7, 4, 5, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7, 5, 9], [7, 4, 5, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([4, 5, 6], [7, 4, 5, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([9], [7, 4, 5, 6, 9])}')
    print(f'{[7, 4, 5, 6, 9]}--{merge_2_list([7], [7, 4, 5, 6, 9])}')

    # print(f'{[3,5,7]}--{merge_2_list([7],[3,5])}')
    # print(f'{[3,5,7]}--{merge_2_list([3,5],[7])}')
    # print(f'{[3,5]}--{merge_2_list([3,5],[3])}')
    # print(f'{[3,5]}--{merge_2_list([3,5],[5])}')
    # print(f'{[3,5]}--{merge_2_list([3],[3,5])}')
    # print(f'{[3,5]}--{merge_2_list([5],[3,5])}')
    # print(f'{[3,5]}--{merge_2_list([3,5],[3,5])}')
    # print(f'{[3,5]}--{merge_2_list([],[3,5])}')
    # print(f'{[3,5]}--{merge_2_list([3,5],[])}')
    # print(f'{[]}--{merge_2_list([],[])}')

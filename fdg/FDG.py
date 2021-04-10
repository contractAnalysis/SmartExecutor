from copy import copy
import numpy as np
import csv

from fdg.funtion_info import Function_info


class FDG():

    def __init__(self,functions_dict):
        """
        :param functions_dict:
        key: [function_name,read_sv_list,write_sv_list,identifier]
         {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], ...}
        """
        self.limit_maximum_depth=5
        self.ftn_to_index={}
        self.index_to_ftn={}
        self.ftn_to_4bytes={}

        self.labels=set()
        self.label_to_index={}
        self.index_to_label={}
        self.num_ftn=0
        self.num_label=0

        for key, value in functions_dict.items():
            #value[0] is function name (not full name)
            key=key.split('f')
            key_index=int(key[1],10)
            self.ftn_to_index[value[0]]=key_index
            self.index_to_ftn[key_index]=value[0]
            self.ftn_to_4bytes[value[0]]=self.ftn_2_4bytes(value[3])
            self.labels=self.labels.union(value[1]+value[2])
        self.labels = sorted(self.labels)
        for idx,item in enumerate(self.labels):
            self.index_to_label[idx]=item
            self.label_to_index[item]=idx

        self.num_ftn=len(self.ftn_to_index)
        self.num_label=len(self.label_to_index)

        self.sv_read =-np.ones((self.num_ftn,self.num_label))
        self.sv_write = -np.ones((self.num_ftn,self.num_label))
        for key,value in functions_dict.items():
            a=-np.ones(self.num_label)
            b = -np.ones(self.num_label)
            if len(value[1])!=0:
                for label in value[1]:
                    lbl_index=self.label_to_index[label]
                    # a[lbl_index]=1
                    a[lbl_index] = lbl_index

            if len(value[2]) != 0:
                for label in value[2]:
                    lbl_index = self.label_to_index[label]
                    # b[lbl_index] = 1
                    b[lbl_index] = lbl_index
            key = key.split('f')
            key_idx=int(key[1],10)
            self.sv_read[key_idx]=a.T
            self.sv_write[key_idx] = b.T

        # build fdg
        self.matrix_fdg=-np.ones((1,self.num_ftn,self.num_ftn))
        # set to self.num_label if functions have empty sv_read list
        # (constructor can reach them without any dependency,self.num_label as index value is not used)
        for j,item in enumerate(self.sv_read):
            if j>0: # do not consider the case of constructor itself
                if np.array_equal(item,-np.ones(self.num_label)):
                    self.matrix_fdg[0, j, 0] = self.num_label

        for i in range(self.limit_maximum_depth): # the maximum depth to 10
            # in the case of constructor and functions having empty sv_read list
            if i==0:
                write_sv_ary=self.sv_write[0]
                # get label_indices that ftn with index 0 write
                lbl_indices_w=write_sv_ary[write_sv_ary>=0]
                for lbl_idx in lbl_indices_w:
                    # get indices of functions that read label with index lbl_idx
                    ftn_indices=self.sv_read[:,int(lbl_idx)]
                    ftn_indices_r=np.where(ftn_indices>=0)[0]
                    if len(ftn_indices_r)>0:
                        for ftn_idx in ftn_indices_r:
                            if ftn_idx!=0:
                                self.matrix_fdg[0,int(ftn_idx),0]=lbl_idx
                continue

            # in the case of other functions with indices >0
            # get function indices that can be reached from the immediate previous step
            ftn_reached=self.ftn_reached_at_a_depth(-1)

            matrix_i = -np.ones((1, self.num_ftn, self.num_ftn))
            for idx_from in ftn_reached:
                write_sv_ary = self.sv_write[int(idx_from)]
                # get label_indices that ftn with index  write
                lbl_indices_w = write_sv_ary[write_sv_ary >= 0]
                for lbl_idx in lbl_indices_w:
                    # get indices of functions that read label with index lbl_idx
                    ftn_indices = self.sv_read[:, int(lbl_idx)]
                    ftn_indices_r = np.where(ftn_indices >= 0)[0]
                    if len(ftn_indices_r) > 0:
                        for idx_to in ftn_indices_r:
                            if idx_to != idx_from:
                                matrix_i [0, int(idx_to), int(idx_from)] = lbl_idx

            # check the conditions to stop
            if np.array_equal(matrix_i,-np.ones((1, self.num_ftn, self.num_ftn))):
                break
            if i>3 and self.is_contained(matrix_i,self.matrix_fdg[-1,:,:]):
                break
            self.matrix_fdg = np.concatenate((self.matrix_fdg, matrix_i), axis=0)
        self.depth_all_ftns_reached=self.depth_all_nodes_reached()
        self.original_matrix_fdg=copy(self.matrix_fdg)


    def is_contained(self,A,B):
        """
        if B contains A, return True, else return False
        :param A:
        :param B:
        :return:
        """
        assert isinstance(A,np.ndarray)
        assert isinstance(B, np.ndarray)
        c=B-A
        re=c[c<0]
        if len(re)>0:return False
        else: return True


    def depth_all_nodes_reached(self):
        # return the depth (starting with 1) that all function is reached
        mark=np.zeros(self.num_ftn)
        mark[0]=1 # ignore function constructor
        for i in range(self.matrix_fdg.shape[0]):
            ftn_reached=self.ftn_reached_at_a_depth(i)
            for idx in ftn_reached:
                mark[idx]=1
            if len(mark[mark==1])==self.num_ftn:
                return i+1

    def ftn_reached_at_a_depth(self,depth)->list:
        re = self.matrix_fdg[depth, :, :]
        ftn_reached = []
        for j in range(re.shape[0]):
            if len(np.where(re[j, :] >= 0)[0]) > 0:
                ftn_reached.append(j)
        return ftn_reached


    def update(self,depth,row,column,value):
        self.matrix_fdg[depth,row,column]=value

    def next_ftn_names(self,depth,ftn_name):
        ftn_names=[]
        d3=self.ftn_to_index[ftn_name]
        ftn_fdg=self.matrix_fdg[depth,:,d3]
        ftn_idx= np.where(ftn_fdg>=0)[0]
        for idx in ftn_idx:
            ftn_names.append(self.index_to_ftn[idx])
        return ftn_names

    def next_ftn_indices(self, depth, ftn_name):
        ftn_names = []
        d3 = self.ftn_to_index[ftn_name]
        ftn_fdg = self.matrix_fdg[depth, :, d3]
        ftn_idx = np.where(ftn_fdg >= 0)[0]

        return ftn_idx



    def ftn_2_4bytes(self,ftn_identifier:str):
        """
           :param ftn_identifier: 6a7301b8
           :return: ['6a','73','01','b8']
           """
        if ftn_identifier is None: return None
        if len(ftn_identifier) != 8: return None
        ftn_bytes = []

        for i in range(0, len(ftn_identifier), 2):
            ftn_bytes.append(int(ftn_identifier[i:i + 2], 16))
        return ftn_bytes

    def write_CSV(self,file_path):
        # file_path='/media/sf___share_vms/fdg_matrix.csv'
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            header = [' ']
            for c in range(self.matrix_fdg.shape[2]):
                header += [self.index_to_ftn[c]]

            for d in range(self.matrix_fdg.shape[0]):
                writer.writerow(header)
                for i in range(self.matrix_fdg.shape[1]):
                    row = [self.index_to_ftn[i]]
                    for j in range(self.matrix_fdg.shape[2]):
                        if self.matrix_fdg[d, i, j] == -1: # means there is no depencency
                            row += [' ']
                        elif self.matrix_fdg[d, i, j]<self.num_label and self.matrix_fdg[d, i, j] >= 0:
                            row += [self.index_to_label[int(self.matrix_fdg[d, i, j])]]
                        elif self.matrix_fdg[d, i, j] == self.num_label:  # means the edge does not have label, directly reachable from initial state
                            row += ['None']
                        elif self.matrix_fdg[d, i, j]==-3: # to indicate the dependency does not exist
                            row+=['false_dependency']
                        elif self.matrix_fdg[d, i, j]==-2:
                            row+=['miss_dependency']
                        else:
                            row+=['ERROR']

                    writer.writerow(row)
                writer.writerow('\n')






if __name__=='__main__':
    # get data needed to draw graph
    ftn_info = Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken')
    functionsDict=ftn_info.functions_dict_slither()
    print(functionsDict)

    fdg=FDG(functionsDict)
    print(fdg.matrix_fdg)

    # fdg = FDG(Function_info('/home/wei/PycharmProjects/Contracts/_wei/HoloToken.sol', 'HoloToken').functions_dict_slither())
    # print(fdg.matrix_fdg)
    # print(fdg.depth_all_ftns_reached)


    # #====================================
    # functions_dict={'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], 'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'], 'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'], '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'], 'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'], 'f6': ['decreaseApproval', [], ['allowed'], '66188463'], 'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'], 'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'], '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'], 'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'], 'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'], 'f0': ['constructor', [], ['owner'], '8afc3605']}
    # fdg=FDG(functions_dict)
    # print(fdg.depth_all_ftns_reached)

    #
    # a=np.array([2,3,4,5])
    # b=a[a>2]
    # print(a)
    # print(b)

    # print(set(fdg.ftn_reached_at_a_depth(1)).difference(set(fdg.ftn_reached_at_a_depth(2))))

    # fdg.write_CSV('/media/sf___share_vms/fdg_matrix.csv')








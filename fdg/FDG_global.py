from copy import copy
import numpy as np


# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=50

# provide them to FDG_pruner for FDG building
global solidity_path
global contract

global depth_all_ftns_reached
depth_all_ftns_reached=2

# save the coverage (from coverage_plugin)
global coverage
coverage=[0,0]

# get instruction indices for each function (from soliditycontract)
global ftns_instr_indices
ftns_instr_indices={}

# save the lists that record which instruction is covered (from coverage_plugin)
global ftns_instr_cov
ftns_instr_cov=[[],[]]

global mapping
mapping=[]

global solc_indices

global method_identifiers
method_identifiers={}

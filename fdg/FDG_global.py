from copy import copy
import numpy as np


# control the number of symbolic transactions issued by LaserEVM
global transaction_count
transaction_count=10

# provide them to FDG_pruner for FDG building
global solidity_path
global contract

global depth_all_ftns_reached
depth_all_ftns_reached=2

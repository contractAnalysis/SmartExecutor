#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import traceback
from typing import Optional, List

from fdg.funtion_info import Function_info
from mythril.laser.ethereum.iprof import InstructionProfiler
from . import MythrilDisassembler
from mythril.support.source_support import Source
from mythril.support.loader import DynLoader
from mythril.support.support_args import args
from mythril.analysis.symbolic import SymExecWrapper
from mythril.analysis.callgraph import generate_graph
from mythril.analysis.traceexplore import get_serializable_statespace
from mythril.analysis.security import fire_lasers, retrieve_callback_issues
from mythril.analysis.report import Report, Issue
from mythril.ethereum.evmcontract import EVMContract
from mythril.laser.smt import SolverStatistics
from mythril.support.start_time import StartTime
from mythril.exceptions import DetectorNotFoundError
from mythril.laser.execution_info import ExecutionInfo
from fdg.FDG import FDG
log = logging.getLogger(__name__)


functions_dict = {'f1': ['transferOwnership', ['owner'], ['owner'], 'f2fde38b'], 'f2': ['transfer', ['balances', 'mintingFinished'], ['balances'], 'a9059cbb'], 'f3': ['transferFrom', ['balances', 'mintingFinished', 'allowed'], ['allowed', 'balances'], '23b872dd'], 'f4': ['approve', ['mintingFinished'], ['allowed'], '095ea7b3'], 'f5': ['increaseApproval', [], ['allowed'], 'd73dd623'], 'f6': ['decreaseApproval', [], ['allowed'], '66188463'], 'f7': ['setMinter', ['owner'], ['minter'], 'fca3b5aa'], 'f8': ['mint', ['balances', 'minter', 'totalSupply'], ['balances', 'totalSupply'], '40c10f19'], 'f9': ['finishMinting', ['minter'], ['mintingFinished'], '7d64bcb4'], 'f10': ['setDestroyer', ['owner'], ['destroyer'], '6a7301b8'], 'f11': ['burn', ['balances', 'destroyer'], ['balances', 'totalSupply'], '42966c68'], 'f0': ['constructor', [], ['owner'], '8afc3605']}

class MythrilAnalyzer:
    """
    The Mythril Analyzer class
    Responsible for the analysis of the smart contracts
    """

    def __init__(
        self,
        disassembler: MythrilDisassembler,
        requires_dynld: bool = False,
        use_onchain_data: bool = True,
        strategy: str = "dfs",
        address: Optional[str] = None,
        max_depth: Optional[int] = None,
        execution_timeout: Optional[int] = None,
        loop_bound: Optional[int] = None,
        create_timeout: Optional[int] = None,
        enable_iprof: bool = False,
        disable_dependency_pruning: bool = False,
        solver_timeout: Optional[int] = None,
        custom_modules_directory: str = "",
        sparse_pruning: bool = False,
        unconstrained_storage: bool = False,
        parallel_solving: bool = False,
        call_depth_limit: int = 3,
        solver_log: Optional[str] = None,
        fdg_flag=False, #@wei

    ):
        """

        :param disassembler: The MythrilDisassembler class
        :param requires_dynld: whether dynamic loading should be done or not
        :param onchain_storage_access: Whether onchain access should be done or not
        """
        # when contract name is given, dissasembly has only one contract
        # when contract name is not given, dissasembly has the number of contracts in the given solidity file
        if fdg_flag:#@wei
            # self.fdg = FDG(functions_dict)
            file_path = disassembler.contracts[0].input_file
            contract_name = disassembler.contracts[0].name
            function_info=Function_info(file_path,contract_name)
            self.fdg=FDG(function_info.functions_dict_slither())
            # print(f'contract={contract}\nfile_path={file_path}')
        else:
            self.fdg =None

        self.eth = disassembler.eth
        self.contracts = disassembler.contracts or []  # type: List[EVMContract]
        self.enable_online_lookup = disassembler.enable_online_lookup
        self.use_onchain_data = use_onchain_data
        self.strategy = strategy
        self.address = address
        self.max_depth = max_depth
        self.execution_timeout = execution_timeout
        self.loop_bound = loop_bound
        self.create_timeout = create_timeout
        self.disable_dependency_pruning = disable_dependency_pruning
        self.custom_modules_directory = custom_modules_directory
        args.sparse_pruning = sparse_pruning
        args.solver_timeout = solver_timeout
        args.parallel_solving = parallel_solving
        args.unconstrained_storage = unconstrained_storage
        args.call_depth_limit = call_depth_limit
        args.iprof = enable_iprof
        args.solver_log = solver_log

    def dump_statespace(self, contract: EVMContract = None) -> str:
        """
        Returns serializable statespace of the contract
        :param contract: The Contract on which the analysis should be done
        :return: The serialized state space
        """
        sym = SymExecWrapper(
            contract or self.contracts[0],
            self.address,
            self.strategy,
            dynloader=DynLoader(self.eth, active=self.use_onchain_data),
            max_depth=self.max_depth,
            execution_timeout=self.execution_timeout,
            create_timeout=self.create_timeout,
            disable_dependency_pruning=self.disable_dependency_pruning,
            run_analysis_modules=False,
            custom_modules_directory=self.custom_modules_directory,
            fdg=self.fdg #@wei here
        )

        return get_serializable_statespace(sym)

    def graph_html(
        self,
        contract: EVMContract = None,
        enable_physics: bool = False,
        phrackify: bool = False,
        transaction_count: Optional[int] = None,
    ) -> str:
        """

        :param contract: The Contract on which the analysis should be done
        :param enable_physics: If true then enables the graph physics simulation
        :param phrackify: If true generates Phrack-style call graph
        :param transaction_count: The amount of transactions to be executed
        :return: The generated graph in html format
        """

        sym = SymExecWrapper(
            contract or self.contracts[0],
            self.address,
            self.strategy,
            dynloader=DynLoader(self.eth, active=self.use_onchain_data),
            max_depth=self.max_depth,
            execution_timeout=self.execution_timeout,
            transaction_count=transaction_count,
            create_timeout=self.create_timeout,
            disable_dependency_pruning=self.disable_dependency_pruning,
            run_analysis_modules=False,
            custom_modules_directory=self.custom_modules_directory,
            fdg=self.fdg #@wei
        )
        return generate_graph(sym, physics=enable_physics, phrackify=phrackify)

    def fire_lasers(
        self,
        modules: Optional[List[str]] = None,
        transaction_count: Optional[int] = None,
    ) -> Report:
        """
        :param modules: The analysis modules which should be executed
        :param transaction_count: The amount of transactions to be executed
        :return: The Report class which contains the all the issues/vulnerabilities
        """
        all_issues = []  # type: List[Issue]
        SolverStatistics().enabled = True
        exceptions = []
        execution_info = None  # type: Optional[List[ExecutionInfo]]
        for contract in self.contracts:
            StartTime()  # Reinitialize start time for new contracts
            try:
                sym = SymExecWrapper(
                    contract,
                    self.address,
                    self.strategy,
                    dynloader=DynLoader(self.eth, active=self.use_onchain_data),
                    max_depth=self.max_depth,
                    execution_timeout=self.execution_timeout,
                    loop_bound=self.loop_bound,
                    create_timeout=self.create_timeout,
                    transaction_count=transaction_count,
                    modules=modules,
                    compulsory_statespace=False,
                    disable_dependency_pruning=self.disable_dependency_pruning,
                    custom_modules_directory=self.custom_modules_directory,
                    fdg=self.fdg #@wei
                )
                issues = fire_lasers(sym, modules)
                execution_info = sym.execution_info
            except DetectorNotFoundError as e:
                # Bubble up
                raise e
            except KeyboardInterrupt:
                log.critical("Keyboard Interrupt")
                issues = retrieve_callback_issues(modules)
            except Exception:
                log.critical(
                    "Exception occurred, aborting analysis. Please report this issue to the Mythril GitHub page.\n"
                    + traceback.format_exc()
                )
                issues = retrieve_callback_issues(modules)
                exceptions.append(traceback.format_exc())
            for issue in issues:
                issue.add_code_info(contract)

            all_issues += issues
            log.info("Solver statistics: \n{}".format(str(SolverStatistics())))

        source_data = Source()
        source_data.get_source_from_contracts_list(self.contracts)

        # Finally, output the results
        report = Report(
            contracts=self.contracts,
            exceptions=exceptions,
            execution_info=execution_info,
        )
        for issue in all_issues:
            report.append_issue(issue)

        return report

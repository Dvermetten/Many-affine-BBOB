import numpy as np
import ioh
from itertools import product
from functools import partial
from multiprocessing import Pool, cpu_count

import sys
import argparse
import warnings
import os

import nevergrad as ng
import nevergrad.common.typing as tp

import pandas as pd

from nevergrad.optimization.optimizerlib import (
    DifferentialEvolution,
    RCobyla,
    DiagonalCMA,
)

from modcma import ModularCMAES
from modde import ModularDE
# from modcma import Parameters, AskTellCMAES

import time

rootname = "/datanaco/vermettendl/Many_affine/"

def runParallelFunction(runFunction, arguments):
    """
        Return the output of runFunction for each set of arguments,
        making use of as much parallelization as possible on this system

        :param runFunction: The function that can be executed in parallel
        :param arguments:   List of tuples, where each tuple are the arguments
                            to pass to the function
        :return:
    """
    

    arguments = list(arguments)
    p = Pool(min(200, len(arguments)))
    results = p.map(runFunction, arguments)
    p.close()
    return results

algs_considered = ['DiagonalCMA', 'modcma', 'DifferentialEvolution','RCobyla', 'modde']

scale_factors = [11. , 17.5, 12.3, 12.6, 11.5, 15.3, 12.1, 15.3, 15.2, 17.4, 13.4,
       20.4, 12.9, 10.4, 12.3, 10.3,  9.8, 10.6, 10. , 14.7, 10.7, 10.8,
        9. , 12.1]

class NG_Evaluator():
    def __init__(self, optimizer):
        self.alg = optimizer
    
    def __call__(self, func):
        np.random.seed(int(time.time()))
        parametrization = ng.p.Array(shape=(func.meta_data.n_variables,)).set_bounds(-5, 5)
        if self.alg == 'modcma':
            c = ModularCMAES(func, d=func.meta_data.n_variables, bound_correction='saturate',
                             budget=int(2000*func.meta_data.n_variables),
                             x0=np.zeros((func.meta_data.n_variables, 1)), target = 1e-12)
            c.run()
        elif self.alg == 'modde':
            c = ModularDE(func, bound_correction='saturate',
                             budget=int(2000*func.meta_data.n_variables),
                             mutation_base = 'target', 
                              mutation_reference= 'pbest',  lpsr=True, 
                              lambda_ = 18*func.meta_data.n_variables, 
                              use_archive= True, 
                              adaptation_method_F= 'shade', 
                              adaptation_method_CR=  'shade')
            c.run()
        else:
            if self.alg in ['ConfiguredPSO', 'EMNA', 'DifferentialEvolution']:
                optimizer = eval(f"{self.alg}")()(
                    parametrization=parametrization, budget=int(2000*func.meta_data.n_variables)
                )
            else:
                optimizer = eval(f"{self.alg}")(
                    parametrization=parametrization, budget=int(2000*func.meta_data.n_variables)
                )
            optimizer.minimize(func)

        
class ManyAffine():
    def __init__(self, weights, instances, opt_loc=1, dim = 5, sf_type = 'min_max'):
        self.weights = weights / np.sum(weights)
        self.fcts = [ioh.get_problem(fid, int(iid), dim) for fid, iid in zip(range(1,25), instances)]
        self.opts = [f.optimum.y for f in self.fcts]
        self.scale_factors = scale_factors
        if type(opt_loc) == int:
            self.opt_x = self.fcts[opt_loc].optimum.x
        else:
            self.opt_x = opt_loc

    def __call__(self, x):
        raw_vals = np.array([ np.clip(f(x+f.optimum.x - self.opt_x)-o,1e-12,1e20) for f, o in zip(self.fcts, self.opts)])
        weighted = (np.log10(raw_vals)+8)/self.scale_factors * self.weights
        return 10**(10*np.sum(weighted)-8)
    
    
def run_ng_optimizer(temp):
    print(temp)
    idx, algname, dim = temp
    
    weights = pd.read_csv("~/affine_bbob/weights.csv", index_col=0)
    iids = pd.read_csv("~/affine_bbob/iids.csv", index_col=0)
    opt_locs = pd.read_csv("~/affine_bbob/opt_locs.csv", index_col=0)
    
    alg = NG_Evaluator(algname)

    f_new = ManyAffine(np.array(weights.iloc[idx]), 
                       np.array(iids.iloc[idx]), 
                       np.array(opt_locs.iloc[idx])[:dim], dim)
    
    ioh.problem.wrap_real_problem(
        f_new,                                     
        name=f"affine_{idx}",                               
        optimization_type=ioh.OptimizationType.MIN, 
        lb=-5,                                               
        ub=5
    )
    
    f = ioh.get_problem(f"affine_{idx}", 0, dim)
    
    l = ioh.logger.Analyzer(root = rootname, 
                            folder_name = f"{idx}_{algname}_{dim}D", 
                            algorithm_name=algname, store_positions=(dim == 2))
    # l.set_experiment_attributes({'weights' : np.array(weights.iloc[idx]), 
                                 # 'iids' : np.array(iids.iloc[idx]), 
                                 # 'opt_loc' : np.array(opt_locs.iloc[idx])[:dim]})
    f.attach_logger(l)
    for rep in range(50):
        np.random.seed(rep)
        alg(f)
        f.reset()
        
    l.close()
    
def run_ng_optimizer_default(temp):
    print(temp)
    fid, iid, algname, dim = temp
    
    f_base = ioh.get_problem(fid, iid, dim)
    
    weights = np.zeros(24)
    weights[fid-1] = 1
    
    iids = np.ones(24)*iid
    
    f_new = ManyAffine(weights, 
                       iids, 
                       f_base.optimum.x, dim)
    
    alg = NG_Evaluator(algname)
    
    ioh.problem.wrap_real_problem(
        f_new,                                     
        name=f"F{fid}_I{iid}",                               
        optimization_type=ioh.OptimizationType.MIN, 
        lb=-5,                                               
        ub=5
    )
    
    f = ioh.get_problem(f"F{fid}_I{iid}", 0, dim)
    
    

    l = ioh.logger.Analyzer(root = rootname, 
                            folder_name = f"F{fid}_I{iid}_{algname}_{dim}D", 
                            algorithm_name=algname, store_positions=(dim == 2))
    # l.set_experiment_attributes({'weights' : np.array(weights.iloc[idx]), 
                                 # 'iids' : np.array(iids.iloc[idx]), 
                                 # 'opt_loc' : np.array(opt_locs.iloc[idx])[:dim]})
    f.attach_logger(l)
    for rep in range(50):
        np.random.seed(rep)
        alg(f)
        f.reset()
        
    l.close()
    

if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    warnings.filterwarnings("ignore", category=FutureWarning)

    idxs = range(1000)
    algnames =  ['modcma', 'DifferentialEvolution', 'DiagonalCMA', 'modde', 'RCobyla']
    dims = [2,5]
    args = product(idxs, algnames, dims)
    
    runParallelFunction(run_ng_optimizer, args)

    fids = range(1,25)
    iids = range(1,6)
    dims = [2,5]
    args = product(fids, iids, algnames, dims)
    
    runParallelFunction(run_ng_optimizer_default, args)
    
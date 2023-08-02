import numpy as np
import ioh
from itertools import product
from functools import partial
from multiprocessing import Pool, cpu_count

import sys
import argparse
import warnings
import os

from scipy.stats import qmc

import pandas as pd

import pflacco.classical_ela_features as pflacco_ela

dirname = "/mnt/d/Many_Affine/ELA2"

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
    p = Pool(min(8, len(arguments)))
    results = p.map(runFunction, arguments)
    p.close()
    return results

scale_factors = [11. , 17.5, 12.3, 12.6, 11.5, 15.3, 12.1, 15.3, 15.2, 17.4, 13.4,
       20.4, 12.9, 10.4, 12.3, 10.3,  9.8, 10.6, 10. , 14.7, 10.7, 10.8,
        9. , 12.1]

def Sobol_sampling(dimension, sample_size, random_seed):
    sampler = qmc.Sobol(d=dimension, scramble=True, seed=random_seed)
    # sample = sampler.random_base2(m=3)
    sample = sampler.random(n=sample_size)
    return sample

def compute_ela(X, y, lower_bound=-5, upper_bound=5):
    y_rescale = (max(y) - y) / (max(y)-min(y))
    # Calculate ELA features
    ela_meta = pflacco_ela.calculate_ela_meta(X, y_rescale)
    ela_distr = pflacco_ela.calculate_ela_distribution(X, y_rescale)
    ela_level = pflacco_ela.calculate_ela_level(X, y_rescale)
    pca = pflacco_ela.calculate_pca(X, y_rescale)
    limo = pflacco_ela.calculate_limo(X, y_rescale, lower_bound, upper_bound)
    nbc = pflacco_ela.calculate_nbc(X, y_rescale)
    disp = pflacco_ela.calculate_dispersion(X, y_rescale)
    ic = pflacco_ela.calculate_information_content(X, y_rescale, seed=100)
    ela_ = {**ela_meta, **ela_distr, **ela_level, **pca, **limo, **nbc, **disp, **ic}
    df_ela = pd.DataFrame([ela_])
    return df_ela
        
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
    
    
def calc_ela_affine(temp):
    print(temp)
    idx, dim = temp
    
    weights = pd.read_csv("weights.csv", index_col=0)
    iids = pd.read_csv("iids.csv", index_col=0)
    opt_locs = pd.read_csv("opt_locs.csv", index_col=0)
    
    f_new = ManyAffine(np.array(weights.iloc[idx]), 
                       np.array(iids.iloc[idx]), 
                       np.array(opt_locs.iloc[idx])[:dim], dim)
    
    df_ela_prob = pd.DataFrame()
    for i in range(5):
        doe_x = Sobol_sampling(dim, dim*1000, i)
        doe_x = doe_x*10-5
        y = np.array(list(map(f_new, doe_x)))
        df_ela = compute_ela(doe_x, y, lower_bound=-5, upper_bound=5)
    # df_ela['ndoe'] = i
        df_ela_prob = pd.concat([df_ela_prob, df_ela], axis=0, ignore_index=True)
        
    df_ela_prob.to_csv(f"{dirname}/{idx}_{dim}D.csv")
    
def calc_ela_default(temp):
    print(temp)
    fid, iid, dim = temp
    
    f_base = ioh.get_problem(fid, iid, dim)
    
    weights = np.zeros(24)
    weights[fid-1] = 1
    
    iids = np.ones(24)*iid
    
    f_new = ManyAffine(weights, 
                       iids, 
                       f_base.optimum.x, dim)
    
    df_ela_prob = pd.DataFrame()
    for i in range(5):
        doe_x = Sobol_sampling(dim, dim*1000, 42)
        doe_x = doe_x*10-5
        y = np.array(list(map(f_new, doe_x)))
        df_ela = compute_ela(doe_x, y, lower_bound=-5, upper_bound=5)
        df_ela['ndoe'] = i
        df_ela_prob = pd.concat([df_ela_prob, df_ela], axis=0, ignore_index=True)
        
    df_ela_prob.to_csv(f"{dirname}/F{fid}_I{iid}_{dim}D.csv")
    
    
if __name__ == '__main__':
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    warnings.filterwarnings("ignore", category=FutureWarning)

    idxs = range(1000)
    dims = [2,5]
    args = product(idxs, dims)
    
    runParallelFunction(calc_ela_affine, args)
    
    fids = range(1,25)
    iids = range(1,6)
    dims = [2,5]
    args = product(fids, iids, dims)
    
    runParallelFunction(calc_ela_default, args)
    
    


    
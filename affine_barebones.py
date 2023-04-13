import ioh
import numpy as np

scale_factors = [11. , 17.5, 12.3, 12.6, 11.5, 15.3, 12.1, 15.3, 15.2, 17.4, 13.4,
       20.4, 12.9, 10.4, 12.3, 10.3,  9.8, 10.6, 10. , 14.7, 10.7, 10.8,
        9. , 12.1]

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
    
    
if __name__ == '__main__':
    #This is an example of how to access the function generator to make a 5-dimensional problem
    
    dim = 5
    
    weights = np.random.uniform(size=24)
    iids = np.random.randint(100, size=24)
    opt_loc = np.random.uniform(size=dim)*10-5 #in [-5,5] as BBOB
    
    f_new = ManyAffine(weights, 
                       iids, 
                       opt_loc, dim)
    
    ioh.problem.wrap_real_problem(
        f_new,                                     
        name=f"F{fid}_I{iid}",                               
        optimization_type=ioh.OptimizationType.MIN, 
        lb=-5,                                               
        ub=5
    )
    
    f = ioh.get_problem(f"F{fid}_I{iid}", 0, dim)
    
    #Now, f is usable as any regular problem from IOHexperimenter, e.g. for adding loggers
    
    print(f.meta_data)

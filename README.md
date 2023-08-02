This document details the reproduction steps for the paper "Many-affine BBOB function generator"

# Dependencies

This project relies on a few different libraries to achieve the presented results. These are as follows:

## IOHexperimenter

Part of the [IOHprofiler](https://iohprofiler.github.io/) environment, [IOHexperimenter](https://iohprofiler.github.io/IOHexp/) provides an interface between optimization algorithms and problems, and adds detailed logging functionality to this pipeline. 
We use the python-version of IOHexperimenter, available on [pip as 'ioh'](https://pypi.org/project/ioh/) (we used version 0.3.9). From this package, we use the logging component, as well as the interface to the [BBOB problem suite](https://bee22.com/resources/bbob%20functions.pdf). 

## PFlacco

To analyze the function's low-level properties, we make use of Exploratory Landscape Analysis (ELA), which gives access to a wide range of features. To calculate these, we use the python-based [pflacco](https://github.com/Reiyan/pflacco) library (version 1.2.0, note that this requires python 3.8 or higher).  
We make use of only the features which don't require sampling additional points from the function, from the 'classical_ela_features' module. 

## Nevergrad

To access a variety of optimization algorihms, we make use of the [Nevergrad](https://github.com/facebookresearch/nevergrad) (version 0.4.3.post8). 
We make use of the following algorithms from Nevergrads optimizers module:  'DifferentialEvolution', 'DiagonalCMA', 'RCobyla'

## Modular CMA-ES + DE

In addition to the Nevergrad algorithms, we make use of two modular algorithms frameworks in our portfolio. The first is [Modular CMA-ES](https://github.com/IOHprofiler/ModularCMAES), 'modcma' on pip (version 0.0.2.8.4). The second is [Modular DE](https://github.com/Dvermetten/ModDE), 'modde' on pip (version 0.0.1).

## IOHanalzyer

As a final requirement, we make use of the [IOHanalyzer](https://github.com/IOHprofiler/IOHanalyzer). This is an R-based library for analyzing and visualizing optimization algorithm performance. We use version 0.1.7.2.

# Core files in this repositoy

The most important functionality of this paper is of course the MA-BBOB function generator. Because of this, even as all experiment scripts have the implementation of the function generator included, we also add a barebones version of the generator in 'affine_barebones.py' to show how to use it as a normal iohprofiler-based benchmark problem.

This 'affine_barebones.py' shows in its main function how we can create a arbitrary instance of MA-BBOB, and use the IOHexperimenter to wrap it into a regular 'ioh' function object. For more information on how to work with this type of object, we refer to the [IOHexperimenter tutorial](https://github.com/IOHprofiler/IOHexperimenter/blob/master/example/tutorial.ipynb), which shows how to use this function object with any optimization algorithm and combine it with the wide variety of logging functionality.

# Reproducing the papers results

The core file for reproducing the results from this project is the notebook 'Visualization.ipynb'. This notebook is interrupted at times to run scripts and collect data, as we will discuss below. 

## Initial exploration: scaling factor and 2D figures

Before running large-scale expleriments, we start by investigating how we can combat the scaling problem. Since each BBOB function can have a widely different scale, we perform random samplings on each function to gain an estimate of the function value ranges. 
Note that we care about the ranges in relation to the global optimum value, so we only look at the 'precision' here. 

With the random samples collected, we add the implementation of the MA-BBOB generator and add code to plot the landscapes in 2D. This is used to expleriment with different ways to set the scale factors, optimum locations and weight schemes. 

## Determine the settings used
The section 'Setup performance data collection' is used to generate the used weights, instance number and optima location. Each of these are stored to their corresponding csv-file.
These files fully specify the 1000 instances of the MA-BBOB suite we use throughout the paper. These files are used in the scripts for the following parts (ELA + Performace) 

## Calculating ELA features
The ELA-feature computation is based on pflacco as described in the dependencies. 
The script: 'ELA_computation_local.py' runs the computation and stores the results as csv files (Note: the 'dirname' parameter should be changed before running this script). The file loops over all selected BBOB and MA-BBOB instances and gets the following sets of ELA features:
*    meta data
*    distribution
*    level set
*    principal component analysi
*    linear model
*    nbc
*    dispersion
*    information_content
For more information on these feature sets, please see "Mersmann et al. (2011), “Exploratory Landscape Analysis”, in Proceedings of the 13th Annual Conference on Genetic and Evolutionary Computation, pp. 829—836. ACM (http://dx.doi.org/10.1145/2001576.2001690)"

## Performance data

### Collect performance data
The performance data collection script is 'data_collection_manyaffine.py', which makes use of nevergrad and iohprofiler to benchmark the selected algorithms. Note that the 'rootname' parameter should be modified when running this script.

The data from running this script is available in the zenodo repository as 'data.zip'.

### Process performance data
The previous script generates IOHanalyzer-compatible data. This can be processed via the R-package 'IOHanalyzer' using the script 'processing.R'. This results in the files in 'csvs.zip' (in the zenodo), of which we use the auc-based ones which are concatenated into 'aucs.csv'.

## Visualization and algorithm selector
All remaining analysis and visualization is part of the notebook. Note that the corresponding directory names should be updated according to the ones used in the respective scripts. 
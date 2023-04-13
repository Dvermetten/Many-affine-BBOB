This document details the reproduction steps for the paper "Many-affine BBOB function generator"

The core file for this project is the notebook 'Visualization.ipynb'

# The affine function generator

While all experiment scripts have the implementation of the function generator included, we also add a barebones version of the generator in 'affine_barebones.py' to show how to use it as a normal iohprofiler-based benchmark problem.

## Determine the settings used
In the notebook, the section 'Setup performance data collection' is used to generate the used weights, instance number and optima location. Each of these are stored to their corresponding csv-file.

## Calculating ELA features
The ELA-feature computation relies on the 'pflacco' package, which requires python>=3.8. The script: 'ELA_computation_local.py' runs the computation and stores the results as csv files (the directory should be changed before running)

## Performance data

### Collect performance data
The performance data collection script is 'data_collection_manyaffine.py', which makes use of nevergrad and iohprofiler to benchmark the selected algorithms. The data is available in the corresponding zenodo repository as 'data.zip'.

### Process performance data
The previous script generates IOHanalyzer-compatible data. This can be processed via the R-package 'IOHanalyzer' using the script 'processing.R' This results in the files in 'csvs.zip' (in the zenodo), of which we use the auc-based ones which are concatenated into 'aucs.csv'.

## Visualization and algorithm selector
All remaining analysis and visualization is part of the notebook
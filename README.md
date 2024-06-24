# Controlled Noise: Evidence for the Epigenetic Control of Single-Cell Expression Variability

## Data Availability

Our data files are too big to upload to github. We upload the preprocessed and homogeneity-checked data to `Releases`. 

Each `Releases` includes data that has been preprocessed and homogeneity checked, a list of index vectors indicating nearby peaks of genes, index used in train-test splits, and data for Python code.

More detailed infomation can be found in `Releases` description.

`HSPC_Donor31800`: 3000 cells, 206611 peaks and 9657 genes.

`HSPC_Donor32606`: 3000 cells, 207844 peaks and 9711 genes.

`Neurons`: 2000 cells, 26330 peaks and 3003 genes.


## Files Description

`functions_used.R`: R functions used in regression. You can use it by `source("functions_used.R")`.

`make_predictions.R`: After importing functions and data, use Lasso, PeakFusion, LocReg, KNN and GS to make predictions.

`LGBM_training.py`: One train-test split of LGBM.

`MNN_training.py`: One train-test split of MNN.

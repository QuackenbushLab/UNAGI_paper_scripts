# UNAGI Paper Replication

## Download Gene Expression Data
Download gene expression data from GRAND for the following data sets, saving all to the same directory:
- [Skeletal muscle](https://granddb.s3.amazonaws.com/tissues/expression/Skeletal_Muscle.csv)
- [Subcutaneous adipose](https://granddb.s3.amazonaws.com/tissues/expression/Adipose_Subcutaneous.csv)
- [Skin](https://granddb.s3.amazonaws.com/tissues/expression/Skin.csv)
- [Lung](https://granddb.s3.amazonaws.com/tissues/expression/Lung.csv)
- [Aorta](https://granddb.s3.amazonaws.com/tissues/expression/Artery_Aorta.csv)

## Install netZooR
You will need to install the GitHub version of netZooR to run UNAGI. Follow the [instructions to install netZooR](https://github.com/netZoo/netZooR).

## Generate Co-Expression Networks
Run `GenerateCoexpression.R`, changing the following variables in the script prior to running:
   - `expressionDir`: Directory where expression data are stored
   - `coexpressionDir`: Directory where you wish to store coexpression data

## Run UNAGI on Co-Expression Networks
Run `UNAGI_VALIDATION.R`, changing the following variables in the script prior to running. This script generates Supplementary Tables 1-5 as well as filtered coexpression data.
   - `coexpressionDir`: Directory where coexpression data are stored
   - `filteredCoexpressionDir`: Directory where you wish to store the filtered coexpression
   - `unagiDir`: Directory where you wish to store UNAGI results

## Generate Randomized Gene Sets
Run `generate_rand_sets.R`, changing the following variables in the script prior to running:
   - `filteredCoexpressionDir`: Directory where the filtered coexpression is stored
   - `randomSetDir`: Directory where you wish to store the random gene sets

## Run UNAGI on Randomized Gene Sets
We run UNAGI on randomized gene sets to compare the connectivity between actual genes of interest to connectivity between random genes. To do this, run `unagi_serial_rand_runs.R`, changing the following variables in the script prior to running:
   - `coexpressionDir`: Directory where coexpression data are stored
   - `randomSetDir`: Directory where the random gene sets are stored
   - `randResultDir`: Directory where the UNAGI results on random gene sets should be stored

## Analyze UNAGI Results
To plot the networks in Figure 1 and Supplementary Figure 1 and to obtain the results in Tables 1 and 2, run `UNAGI_GRAND_analysis.R`, changing the following variables in the script prior to running:
   - `unagiDir`: Directory where UNAGI results are stored
   - `randResultDir`: Directory where the UNAGI results on random gene sets are stored
   - `tableDir`: Directory where result tables should be stored
   - `plotDir:` Directory where network plots should be stored



# A recursive binning algorithm to measure statistical dependence

The code has been split into three files:
- **binMethods.R** provides the formal definition of a bin object, which is used for convenience in processing bins and flattening the recursive bin structure
- **recursiveBins.R** has the central functions of the algorithm, including the scoring, splitting, and recursive wrapper functions
- **maximizedBinning.R** contains the experimentation used to derive the null distribution, investigate different parameters, and explore the abalone data set

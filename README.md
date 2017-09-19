# Breast cancer CYCLOPS

generate CYCLOPS genes for breast cancer
procedure from Nijhawan et al. (2012)

## Copy Number Analysis of Cancer Cell Lines from CCLE

1. For each sample, creat 100 bin histogram of copy number values for all markers.
2. Use 5-bin moving average to smooth distribution.
- This typically yields 2-5 well-seperated peaks (local maxima with height as measured from local max to rurrounding local minima >2% of genome), presumably corresponding to integer level copy loss and gains.
3. Based on these peaks, seperate samples into one of two categories for classficiation.
- If one peak between log2 copy number -0.05 and -.-5, with as second peak between -0.05 and 0.4 -> First peak copy neutral, second peak partial copy loss. In this case, the cutoff for copy loss was set at 95% upper bound of the second peak, and the cutoff for copy neutral was set at 95% lower bound of the first peak.
- If no peak that meets the criteria, markers <-0.4 -> copy loss, markers >-0.2 -> copy neutral.
4. For markers that lay between the two cutoffs were left uncalled and removed from further anaysis.
5. Marks with log2 copy number ratios <= -1.28 were considered homozygous loss and genew with these copy numbers were also removed.

## Lab Note by date

#### 2017-09-18

* organized ipynb file to python source code
* implemented peak finding

##### 2017-09-19

* working on determining 95% bounds
* working on determining copy number neutral / loss
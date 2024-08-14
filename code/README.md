This directory contains all code required to obtain the results shown Tables 1, 2, 3, and 4, in the main article. 

data-gen-simulation-code.R  : used to generate data as per simulation plan described in Section 5

alasso-simulation-code.R  : used to obtain ALASSO based boostrap confidence intervals by using the simulated data. Relevant results are shown in Table 2.

post-lasso-ols-simulation-code.R  : used to obtain Post-Lasso OLS based boostrap confidence intervals by using the simulated data. Relevant results are shown in Table 1.

combined-real-data-analysis.R  :  does reading, preprocessing and computation relevant confidence intervals from the real data set Ro131.csv. Relevant results are shown in Tables 3 and 4. This code internally calls alasso-real-data.R and post-lasso-ols-real-data.R.

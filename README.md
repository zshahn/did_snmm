# did_snmm
Code and data for applications from Structural Nested Mean Models Under Parallel Trends Assumptions

Apologies, this code is messy! I will work on making this neater over time. I'm also working with a much better programmer to make an R package to implement these methods.

did_snmm_bank is a zipped folder containing R code and data used for the bank deregulation analyses in the paper. It was previously analyzed by Favara and Imbs (2015) and Chaisemartin and D'haultfoueille (2021). 

did_snmm_flood contains code for analyzing the flood dataset from Gallagher (2015). The data itself is too large to fit on github but can be obtained from IPCSR: https://www.openicpsr.org/openicpsr/project/113898/version/V1/view

did_snmm_reboot.pdf is an updated draft of the manuscript that has not yet been posted on arxiv. It contains some additional details (see Remarks 2 and 6) on implementing the machine learning estimation approach described in the paper. Neither real data application employed machine learning because there were not enough occurrences of treatment at each time point. 

The file did_snmm_crossfitting.R contains code for a simulation implementing machine learning nuisance model estimation with cross-fitting as described in Section 3.2 of the paper through Remark 2.

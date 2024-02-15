# GSR_PID
code for simulation and main analysis

This code is used in the paper "Disambiguating the role of blood flow and global signal with partial information decomposition".
https://doi.org/10.1016/j.neuroimage.2020.116699

There are two main files:

1. running the PID oh HCP data, with 278 ROIs, the Global Signal, and the BOLD from a vessel mask.

The PID approach is described here "Multiscale Analysis of Information Dynamics for Linear Multivariate Processes",
by Faes, Marinazzo and Stramaglia http://www.mdpi.com/1099-4300/19/8/408, 2017 Entropy, and the code is here https://github.com/danielemarinazzo/multiscale_PID

2. simulating "neural" and "vessel" BOLD time series with simTB http://mialab.mrn.org/software/simtb/documentation.html

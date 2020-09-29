These are the R scripts accompanying Bartoszek, Vasterlund "'Old Techniques for New Times': the RMaCzek package for producing Czekanowski’s Diagrams".

The R setup for the manuscript was as follows:
R version 3.6.1 (2019-09-12)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: openSUSE Leap 42.3

The exact output can depend on the random seed.

The code is divided into several scripts.
1) simulations_Section3.R
    Code for the clustered phylogeny simulation study of Section 3.
    The Random_seed_1.RData and Random_seed_plot.RData files contain the random seeds to replicate these results.
    This script can take a long time to run, on the order of 4 days (3.50GHz Intel(R) Xeon(R) CPU). To reduce
    the running time one can decrease the value of the numreps variable or the entries in the v_clustsize
    vector (defined after function definitions).

2) JCze1909.R
    Code for the reanalyses of Czekanowski (1909)'s skull data in Section 5. 
    The exact output can depend on the random seed due to possibly different trajectories of the seriation methods.
    The file JCze1909_randomseed.RData contains the random seed to replicate the results.

3) ASolPJas1999.R 
    Code for the reanalyses of Soltysiak and Jaskulski (1999)'s urns data in Section 5.
    The exact output can depend on the random seed due to possibly different trajectories of the seriation methods.
    The file ASolPJas1999_randomseed.RData contains the random seed to replicate the results.

4) ASol2000.R 
    Code for the reanalyses of Soltysiak (2000)'s seals data in Section 5.
    The exact output can depend on the random seed due to possibly different trajectories of the seriation methods.
    The file ASol2000_randomseed.RData contains the random seed to replicate the results.

5) KWar2015.R 
    Code for the reanalyses of Warzecha (2015)'s internet availability data in Section 5.
    The exact output can depend on the random seed due to possibly different trajectories of the seriation methods.
    The file KWar2015_randomseed.RData contains the random seed to replicate the results.


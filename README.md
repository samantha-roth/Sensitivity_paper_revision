To run the analysis, simply follow the scripts based on their beginning numbers:

0_library.R is the script that loads all the packages, functions and parameters. You may change the test function to use or the testing dimensions or the sample sizes in this script. (you also define the working directory here, essentially you only to make changes in this script)

1_ through 4_ scripts are the analysis of the four sensitivity analysis methods. They record the convergence size and sensitivity indices under different sample sizes. 

5_ scripts prepare the data needed for figures 3 and 4; 6_ scripts prepare the data needed for figures 5 and 6. You may run both 5_ scripts after finishing 1_ through 4_, similarly after running both 5_ scripts, you can run both 6_ scripts.

7_ scripts generate the figures 3 through 6.

8_Gridplot.R generates figures 7 and 8.

9_Supplementary.R generates the supplementary figures.

Hymod.R and sacSma.R are the two hydro models. arnosubbiano.rda and SacSma_dataframe are the corresponding data required to run the two hydro models. sobol_indices_boot.R is a function that allows bootstrapping using sensobol package.

The folder Ranking_Data includes all the required data. You may compare your running results with these data, but notice that the recorded sample sizes, running time and sensitivity indices won't be the same because of random seeds and different computational environments. However you should get similar figures (compare with the figures in New_Figures folder).

**Note that 2_, 3_ and 4_ scripts take a very long time to run, you may begin with low-dimensional test functions to make sure the script can run successfully.

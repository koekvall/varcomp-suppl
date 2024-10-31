# varcomp-suppl
Supplementary material for ["Fast and reliable confidence intervals for a
variance component"](https://arxiv.org/pdf/2404.15060) by Zhang, Ekvall,
and Molstad.

# Figure 1 (Examples)
To reproduce Figure 1 from the main manuscript and Figure 1 from the
Supplementary Material, execute the CIsExample_R1.R and CIsExample_R1_Zoomed.R
scripts from the "Examples" directory. The corresponding .sh files (SLURM job
submission scripts written in Bash) can be used to run these on a
high-performance computing cluster. Once the jobs are complete, the plots can be
created by executing CoveragePlots_R1.R. 

# Figure 2 (Convexity)
Figure 2 can be created by executing the script quasiconvex.R in the "Convexity"
directory. 

# Figure 3 (EigenSim)
The results displayed in Figure 3 can be reproduced by running the script
EigenSim_R1.R for all combinations of "params" (line 5). This can be done in
parallel by executing "EigenSim_R1.sh" (a SLURM job submission script written in
Bash) on a high-performance computing cluster. Once the results are available
(and stored in a directory Results_R1), the plots can be created by executing
the EigenSimPlots_R1.R script. 


# Figure 4  (TimingComparison)
The results displayed in Figure 4 can be produced by using the files available
in the directory "TimingComparison". Specifically, one must first run
TimingSetup_R1.R and store results in a directory "data_R1". Then, one must run
the scripts TimingComparison1_R1.R, TimingComparison01_R1.R, and
TimingComparison5_R1.R one hundred times each, with each of the 100 runs using a
different "uu" (line 5). This can be done automatically by executing the three
corresponding .sh scripts (SLURM job submission scripts written in Bash) on a
high-performance computing cluster. Once the results are available (and stored
in a directory "Results_R1"), the plots can be created by first running the
getResults_1.R script, then running the EigenSimPlots_R1.R script. 

# Figure 5 (SpatialTranscriptomics)
The plots from Figure 5 can be reproduced by running both stxbrain_posterior.R
and stxbrain_anterior.R from the "SpatialTranscriptomics" directory. 

# Figure 2 of Supplementary Materials (WidthComparison)
The results displayed in Figure 2 of the Supplementary Material can be
reproduced by running Widths_R1.R for all combinations of parameters defined in
"params" (line 5). This can be done automatically on a high-performance
computing cluster using the .sh script (SLURM job submission scripts written in
Bash). 


MixedGraphicalModels
====================
Updated: March 28 2014.

This repository contains code to reproduce results in 

Chen, Shizhe, Daniela Witten, and Ali Shojaie. "Selection and Estimation for Mixed Graphical Models." arXiv preprint arXiv:1311.0085 (2013).

Link:
http://arxiv.org/abs/1311.0085

There are three numerical experiments in the paper

1) The empirical probability of successful recovery (Figure 2) in Section 6.2: (./Figure2_prob_recovery).

Code for reproducing Figure 2 can be found in Figure2.r. Note that this code is computationally intensive and might take weeks to run. A toy example can be found in Figure2_quick.r, where we tune down the number of replicates and parameters of the Gibbs sampler. 

2) Edge selection for Gaussian-binary graphs (Figure 3) in Section 6.3: (./Figure3_Gaussian_binary).

Code for reproducing Figure 3 can be found in Figure3.r. Again, A toy example can be found in Figure3_quick.r. Note that the code does not include methods whose code is not publicly available.

3) Application of selection rule on Poisson-binary graphs (Figure 4) in Section 6.4 (./Figure4_Poisson_binary).

Code for reproducing Figure 4 can be found in Figure4.r. Again, A toy example can be found in Figure3_quick.r. Note that the code does not include GRaFo by Fellinghauer et al. 2013.


Functions used in the simulation can be found in ./Sources.

----------
Links to the authors

Shizhe Chen http://students.washington.edu/szchen/
Daniela Witten http://www.biostat.washington.edu/~dwitten/
Ali Shojaie http://www.biostat.washington.edu/~ashojaie/





% File: main.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 06th March 2024
% License: MIT
% Reference:
%{
@article{leeman2024fast,
  title={Fast System Level Synthesis: Robust Model Predictive Control using Riccati Recursions},
  author={Leeman, Antoine P and K{\"o}hler, Johannes and Messerer, Florian and Lahr, Amon and Diehl, Moritz and Zeilinger, Melanie N},
  journal={arXiv preprint arXiv:2401.13762},
  year={2024}}
%}
% Link: https://arxiv.org/abs/2401.13762
% -----------------------------------------------------------------------------
%%
expe00 % simulation with double integrator
expe01 % solvers comparison for increasing horizon length
expe02 % solvers comparison for increasing state dimension
expe03 % evaluation of the number of iteration required before convergence
expe04 % sanity check: all solvers should return the same optimal solution

plots_paper % plots as in the paper
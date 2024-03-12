% File: main.m
% Author: Antoine Leeman (aleeman(at)ethz(dot)ch)
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
addpath('util');

gurobi_installed = true;
mosek_installed = true;
% gurobi and mosek may fail to solve the largest instanciations of the
% problem below depending on the computer used

%% Casadi required for this section ! %
expe00_fast_SLS % simulation with double integrator
expe03_fast_SLS % evaluation of the number of iteration required before convergence

%% solvers comparison for increasing horizon length
expe01_horizon_fast_SLS %casadi required!
if gurobi_installed 
    expe01_horizon_gurobi 
end
if mosek_installed
    expe01_horizon_mosek
end

%% solvers comparison for increasing state dimension
expe02_state_fast_SLS %casadi required!
if gurobi_installed
    expe02_state_gurobi
end
if mosek_installed
    expe02_state_mosek
end

%% sanity check: all solvers should return the same optimal solution
if mosek_installed
    expe04
end

%% Plot all the results

if ~mosek_installed
    disp('Warning: As Mosek is not installed, the related data have not been updated in the plots');
end
if ~gurobi_installed
    disp('Warning: As Gurobi is not installed, the related data have not been updated in the plots');
end
plots_paper % plots as in the paper
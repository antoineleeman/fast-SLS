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
addpath('util/'); % contains utility functions
addpath('sys/'); % contains the system dynamics
addpath('solver/'); % contains the Riccati-based SLS solver, and the yalmip interface to Mosek and Gurobi
addpath('expe/'); % contains the experiements used in the paper
addpath('img/'); % contains the image generated from the experiments, used in the paper
addpath('data/'); % contains the data created by the experiments


gurobi_installed = true;
mosek_installed = true;
casadi_installed = true;
% gurobi and mosek may fail to solve the largest instanciations of the
% problems below depending on the computer used

%% Casadi required for this section ! %
% expe00_fast_SLS % simulation with double integrator (not used in paper)
if casadi_installed
    expe03_fast_SLS % evaluation of the number of iteration required until convergence
end
%% solvers comparison for increasing horizon length
if casadi_installed
    expe01_horizon_fast_SLS
end
if gurobi_installed
    expe01_horizon_gurobi % not tested yet
end
if mosek_installed
    expe01_horizon_mosek
end

%% solvers comparison for increasing state dimension
if casadi_installed
    expe02_state_fast_SLS
end
if gurobi_installed
    expe02_state_gurobi
end
if mosek_installed
    expe02_state_mosek
end

%% sanity check: all solvers should return the same optimal solution
if mosek_installed && casadi_installed
    expe04
end

%% Plot all the results
if ~casadi_installed
    disp('Warning: The proposed method relies on Casadi. Install it to be able to plot new results.');
end
if ~mosek_installed
    disp('Warning: As Mosek is not installed, the related data have not been updated in the plots');
end
if ~gurobi_installed
    disp('Warning: As Gurobi is not installed, the related data have not been updated in the plots');
end

plots_paper % plots as in the paper

% File: expe05.m
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
clear all;
close all;
clc;

L = 4;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
N=10;



kkt = KKT_SLS(N,Q,R,msd,Qf);
solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'mosek');

x0 =4*rand(msd.nx,1)-2;
[feasible,it] = kkt.solve(x0);
feasible = solver_yalmip.solve(x0);


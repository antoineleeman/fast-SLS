% File: expe04.m
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
%
L = 10;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 5;
N=10;

kkt = KKT_SLS(N,Q,R,msd,Qf);
solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'mosek');

it_kkt = [];
faulty_implementation = false;

for ii =1:n_sample
    ii
    x0 =2*rand(msd.nx,1)-1;
    [~,~,~,~,~,V0_kkt] = kkt.solve(x0);
    [~,V0_mosek] = solver_yalmip.solve(x0);

    if norm(full(V0_kkt) - V0_mosek,'inf')>= 1e-4
        norm(full(V0_kkt) - V0_mosek,'inf')
        faulty_implementation = true;
    end
end

if faulty_implementation
    disp('Faulty implementation');
else
    disp('Correct implementation');
end
% File: expe05_rocket2D_fast_SLS.m
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

m = Rocket_2D();
Q = diag([10, 10, 10, 10, 200, 200]);
R = diag([1, 10]);
Qf = Q;
N = 20;

x0 = [10; 0; 0; 0; 12.5; 0];
kkt = KKT_SLS(N,Q,R,m,Qf);
[feasible,it] = kkt.solve(x0);
% File: expe00_fast_SLS.m
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
m = Integrator();
Q = eye(m.nx);
R = eye(m.nu);
Qf = Q;

grid_density = 10;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
N = 15;

kkt = KKT_SLS(N,Q,R,m,Qf); 
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end
%%
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

grid_density = 15;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
%
kkt = KKT_SLS(15,Q,R,m,Qf);
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end






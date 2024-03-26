% File: expe03_fast_SLS.m
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
L = 25;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 1000;
N=25;
X0 =2*rand(msd.nx,n_sample)-1;

solver_kkt = KKT_SLS(N,Q,R,msd,Qf);
it_kkt = [];

for ii =1:n_sample
    ii
    x0 =X0(:,ii);
    tic
    [feasible,it] = solver_kkt.solve(x0);
    time =toc;
    if feasible
        [it_kkt] = [it_kkt;it];
    end
end
disp('percentage solved');
length(it_kkt)/n_sample

save(getUniqueName('it_kkt'),'it_kkt','msd');
save('it_kkt.mat','it_kkt','msd');
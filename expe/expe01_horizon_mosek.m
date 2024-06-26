% File: expe01_horizon_mosek.m
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
init
timings_N_mosek = [];
n_sample = 3;
for nn=1:2:10
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf,'mosek');
    timing_yal = [];
    for ii =1:n_sample
        x0 =X0(:,ii);
        [feasible, ~, time] = solver_yalmip.solve(x0);
            if feasible
                timing_yal = [timing_yal;time];
            end
    end
    if time> 60
        break
    end
    timings_N_mosek = [timings_N_mosek,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_mosek'),'timings_N_mosek','msd')
save('../data/mosek-sls-N.mat','timings_N_mosek','msd')

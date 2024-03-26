% File: expe01_horizon_fast_SLS.m
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
timings_N_exact_kkt = [];
timings_N_exact_ff = [];

IT = [];
for nn=1:10:120
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf);
    timing_kkt = [];
    timing_ricc = [];
    for ii =1:n_sample
        x0 =X0(:,ii);
        [feasible, it, time1,time2] = kkt.solve(x0);
        if feasible
            [timing_kkt] = [timing_kkt;time2+time1];
            timing_ricc = [timing_ricc; time1];
            IT = [IT,it];
        end
    end
    timings_N_exact_kkt = [timings_N_exact_kkt,[nn; mean(timing_kkt);std(timing_kkt)]];
    timings_N_exact_ff = [timings_N_exact_ff,[nn; mean(timing_ricc);std(timing_ricc)]];
end

save(getUniqueName('timings_N_exact_kkt'),'timings_N_exact_kkt','timings_N_exact_ff','msd')
save('data/fast-sls-N.mat','timings_N_exact_kkt','timings_N_exact_ff','msd');
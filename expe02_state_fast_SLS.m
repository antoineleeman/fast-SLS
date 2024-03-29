% File: expe02_state_fast_SLS.m
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
timings_M_kkt = [];
timings_M_kkt_ff = [];

IT_M_kkt = [];
for mm = 1:10:80
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    kkt = KKT_SLS(N,Q,R,msd,Qf);
    timing_mm_kkt = [];
        timing_ff = [];

    for ii=1:n_sample
        x0 =2*rand(msd.nx,1)-1;
        [feasible,it, time1,time2] = kkt.solve(x0);
        if feasible
            timing_mm_kkt = [timing_mm_kkt;time2+time1];
            timing_ff = [timing_ff; time1];
            IT_M_kkt = [IT_M_kkt,it];
        end
    end
    timings_M_kkt = [timings_M_kkt,[msd.nx; mean(timing_mm_kkt);std(timing_mm_kkt)]];
    timings_M_kkt_ff = [timings_M_kkt_ff,[msd.nx; mean(timing_ff);std(timing_ff)]];
end

save(getUniqueName('timings_M_kkt'),'timings_M_kkt','timings_M_kkt_ff','IT_M_kkt','msd','N','n_sample')
save('data/timings_M_kkt.mat','timings_M_kkt','timings_M_kkt_ff','msd','N','n_sample');
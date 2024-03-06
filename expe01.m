% File: expe01.m
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
%%
init
timings_N_exact_kkt = [];
timings_N_exact_ff = [];

IT = [];
for nn=3:10:120
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf);
    timing_kkt = [];
    timing_ricc = [];
    for ii =1:n_sample
        x0 =4*rand(msd.nx,1)-2;
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
save('fast-sls-N.mat','timings_N_exact_kkt','timings_N_exact_ff','msd');

%%
init;
n_sample = 3;
timings_N_gurobi = [];
for nn=3:1:10
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf,'gurobi');
    timing_yal = [];
    for ii =1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
            if feasible
                timing_yal = [timing_yal;time];
            end
    end
        if time> 60
        break
    end
    timings_N_gurobi = [timings_N_gurobi,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_gurobi'),'timings_N_gurobi','msd')
save('gurobi-sls-N.mat','timings_N_gurobi','msd')

%%
init
timings_N_mosek = [];
n_sample = 3;
for nn=3:2:15
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf,'mosek');
    timing_yal = [];
    for ii =1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
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
save('mosek-sls-N.mat','timings_N_mosek','msd')

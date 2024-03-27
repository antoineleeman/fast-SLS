% File: expe02_state_gurobi.m
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
n_sample = 3;
timings_M_gurobi = [];
for mm = 1:1:7
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'gurobi'); 
    timings_mm_yalmip = [];
    for ii=1:n_sample
        x0 =2*rand(msd.nx,1)-1;
        [feasible, ~, time] = solver_yalmip.solve(x0);

        if feasible
            timings_mm_yalmip = [timings_mm_yalmip;time];
        end
    end
   timings_M_gurobi = [timings_M_gurobi,[msd.nx; mean(timings_mm_yalmip);std(timings_mm_yalmip)]];
end

save(getUniqueName('timings_M_gurobi'),'timings_M_gurobi','msd','N','n_sample')
save('timings_M_gurobi.mat','timings_M_gurobi','msd','N','n_sample');
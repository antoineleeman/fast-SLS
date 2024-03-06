% File: expe02.m
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
timings_M_kkt = [];
timings_M_kkt_ff = [];

IT_M_kkt = [];
for mm = 3:10:80
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    kkt = KKT_SLS(N,Q,R,msd,Qf);
    timing_mm_kkt = [];
        timing_ff = [];

    for ii=1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        [feasible,it, time1,time2] = kkt.solve(x0);
        if feasible
            timing_mm_kkt = [timing_mm_kkt;time2+time1];
            timing_ff = [timing_ff; time1];
            IT_M_kkt = [IT_M_kkt,it];
        end
    end
    timings_M_kkt = [timings_M_kkt,[mm; mean(timing_mm_kkt);std(timing_mm_kkt)]];
    timings_M_kkt_ff = [timings_M_kkt_ff,[mm; mean(timing_ff);std(timing_ff)]];
end

histogram(IT_M_kkt)
save(getUniqueName('timings_M_kkt'),'timings_M_kkt','timings_M_kkt_ff','IT_M_kkt','msd','N','n_sample')
save('timings_M_kkt.mat','timings_M_kkt','timings_M_kkt_ff','msd','N','n_sample');


%%
init
n_sample = 3;
timings_M_gurobi = [];
for mm = 3:1:7
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'gurobi'); 
    timings_mm_yalmip = [];
    for ii=1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
        if feasible
            timings_mm_yalmip = [timings_mm_yalmip;time];
        end
    end
   timings_M_gurobi = [timings_M_gurobi,[mm; mean(timings_mm_yalmip);std(timings_mm_yalmip)]];
end

save(getUniqueName('timings_M_gurobi'),'timings_M_gurobi','msd','N','n_sample')
save('timings_M_gurobi.mat','timings_M_gurobi','msd','N','n_sample');


%%
init
timings_M_mosek = [];
n_sample = 3;
for mm = 3:1:7
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'mosek'); 
    timings_mm_yalmip = [];
    for ii=1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
        if feasible
            timings_mm_yalmip = [timings_mm_yalmip;time];
        end
    end
   timings_M_mosek= [timings_M_mosek,[mm; mean(timings_mm_yalmip);std(timings_mm_yalmip)]];
end

save(getUniqueName('timings_M_mosek'),'timings_M_mosek','msd','N','n_sample')
save('timings_M_mosek.mat','timings_M_mosek','msd','N','n_sample');
%%
clear all;
close all;
figure(2);
clf;
load('timings_M_mosek');
load('timings_M_gurobi');
load('timings_M_kkt');

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];

errorbar(timings_M_kkt(1,:), timings_M_kkt(2,:), timings_M_kkt(3,:),'LineWidth',2,'Color', colors(1,:));
hold on;
plot(timings_M_gurobi(1,:), timings_M_gurobi(2,:),'LineWidth',2,'Color', colors(3,:),'Marker','+');
plot(timings_M_mosek(1,:), timings_M_mosek(2,:),'LineWidth',2,'Color', colors(4,:),'Marker','+');
plot(timings_M_kkt(1,:), timings_M_kkt(1,:).^3/500000, 'LineWidth', 2, 'Linestyle',':', 'DisplayName', '$\mathcal{O}(n_x^3)$', 'Color', [.5 .5 .5]);
h = plot(timings_M_mosek(1,:), timings_M_mosek(1,:).^3/50, 'LineWidth', 2, 'Linestyle',':', 'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
h = plot(timings_M_gurobi(1,:), timings_M_gurobi(1,:).^3/10, 'LineWidth',2,'Linestyle',':',  'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('fast-SLS','gurobi','mosek','$\mathcal{O}(n_x^3)$','interpreter','latex');
xlabel('Number of states $n_x$','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;

set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
exportgraphics(gcf,strcat('img/fig2.pdf'),'ContentType','vector');
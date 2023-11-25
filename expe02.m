
%%
clear all;
close all;
clc;

N = 15;
n_sample = 1;

%%
timings_M_kkt = [];
IT_M_kkt = [];
for mm = 3:10:100
    mm
    msd = ChainOfMassSpringDampers(mm);
    Q = eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    kkt = KKT_SLS(N,Q,R,msd,Qf);
    timing_mm_kkt = [];
    for ii=1:n_sample
        x0 =rand(msd.nx,1);
        tic
        [feasible,it] = kkt.solve(x0);
        time =toc;
        if feasible
            timings_mm_kkt = [timing_mm_kkt;time];
            IT_M_kkt = [IT_M_kkt,it];
        end
    end
    timings_M_kkt = [timings_M_kkt,[mm; mean(timings_mm_kkt);std(timings_mm_kkt)]];
end
save(getUniqueName('timings_M_kkt'),'timings_M_kkt','IT_M_kkt','msd','N','n_sample')


%%
timings_M_yal = [];
for mm = 3:3:15
    mm
    msd = ChainOfMassSpringDampers(mm);
    Q = eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf); 
    timings_mm_yalmip = [];
    for ii=1:n_sample
        x0 =rand(msd.nx,1);
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
        if feasible
            timings_mm_yalmip = [timings_mm_yalmip;time];
        end
    end
   timings_M_yal = [timings_M_yal,[mm; mean(timings_mm_yalmip);std(timings_mm_yalmip)]];
end

save(getUniqueName('timings_M_yal'),'timings_M_yal','msd','N','n_sample')
%%
clear all;
close all;
figure(2);
clf;
%load('yalm_scalability_nx.mat');
%load('kkt_scalability_nx.mat');
load('24-Nov-2023_16_30_08__timings_M_yal.mat');
load('24-Nov-2023_16_32_27__timings_M_kkt.mat');

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];


errorbar(timings_M_kkt(1,:), timings_M_kkt(2,:), timings_M_kkt(3,:),'LineWidth',2,'Color', colors(1,:));
hold on;
errorbar(timings_M_yal(1,:), timings_M_yal(2,:), timings_M_yal(3,:),'LineWidth',2,'Color', colors(3,:));

% First plot with DisplayName for legend
plot(timings_M_kkt(1,:), timings_M_kkt(1,:).^3/100000, 'LineWidth', 2, 'DisplayName', '$\mathcal{O}(n_x^3)$', 'Color', [.5 .5 .5]);

% Second plot without legend entry
h = plot(timings_M_yal(1,:), timings_M_yal(1,:).^3/100, 'LineWidth', 2, 'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('fast-SLS','Gurobi','$\mathcal{O}(n_x^3)$','interpreter','latex');
xlabel('Number of states n_x');
ylabel('Computation time [s]');
grid on;
nx_min = [min([timings_M_yal(1,:),timings_M_kkt(1,:)])];
nx_max = [max([timings_M_yal(1,:),timings_M_kkt(1,:)])];
cpu_min = [min([timings_M_yal(2,:),timings_M_kkt(2,:)])];
cpu_max = [max([timings_M_yal(2,:),timings_M_kkt(2,:)])];
axis([nx_min, nx_max, cpu_min/2, cpu_max*2])


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
%exportgraphics(gcf,strcat('img/fig2.pdf'),'ContentType','vector');
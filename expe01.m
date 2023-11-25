%%
clear all;
close all;
clc;

L = 6;
msd = ChainOfMassSpringDampers(L);
Q = eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
grid_density = 5;
% x1_range = linspace(-5,5,grid_density);
% x2_range = linspace(-5,5,grid_density);

timings_N_exact_kkt = [];
%for nn=5:15:100
n_sample = 1;
IT = [];
%%
for nn=3:5:150
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf); %x0 seems unused %% check if the bo are well reset
    timing_kkt = [];

    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
        [feasible,it] = kkt.solve(x0);
        time =toc;
        if feasible
            [timing_kkt] = [timing_kkt;time];
            IT = [IT;it];
        end
    end
    if time> 60
        break
    end
    %     end
    % end
    timings_N_exact_kkt = [timings_N_exact_kkt,[nn; mean(timing_kkt);std(timing_kkt)]];
    
end

histogram(IT)
save(getUniqueName('timings_N_exact_kkt'),'timings_N_exact_kkt','msd')

%%

timings_N_gurobi = [];
for nn=3:3:30
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf,'gurobi');
    timing_yal = [];
    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
            if feasible
                timing_yal = [timing_yal;time];
            end
    end
    timings_N_gurobi = [timings_N_gurobi,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_gurobi'),'timings_N_gurobi','msd')

%%

timings_N_mosek = [];
for nn=5:3:30
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf,'mosek'); %x0 seems unused %% check if the bo are well reset
    timing_yal = [];
    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
        feasible = solver_yalmip.solve(x0);
        time =toc;
            if feasible
                timing_yal = [timing_yal;time];
            end
    end
    if time> 180
        break
    end
    timings_N_mosek = [timings_N_mosek,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_mosek'),'timings_N_mosek','msd')


%%
clear all;
close all;
clf;
%load('24-Nov-2023_15_38_10__timings_N_yal.mat');
load('data/24-Nov-2023_16_45_18__timings_N_yal.mat');
msd.nx

load('24-Nov-2023_21_48_08__timings_N_mosek.mat');
msd.nx
%load('24-Nov-2023_15_36_44__timings_N_exact_kkt.mat');
load('data/24-Nov-2023_16_34_54__timings_N_exact_kkt.mat');
msd.nx

clf;

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];

errorbar(timings_N_exact_kkt(1,:), timings_N_exact_kkt(2,:), timings_N_exact_kkt(3,:),'LineWidth',2,'Color', colors(1,:));
hold on;
errorbar(timings_N_yal(1,:), timings_N_yal(2,:), timings_N_yal(3,:),'LineWidth',2,'Color', colors(3,:));
errorbar(timings_N_mosek(1,:), timings_N_mosek(2,:), timings_N_mosek(3,:),'LineWidth',2,'Color', colors(4,:));

%plot(timings_N_kkt(1,:), timings_N_kkt(1,:).^2/100 ,'LineWidth',2, 'Color', [.5 .5 .5]);
plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(2.5)/1000 ,'LineWidth',2, 'Color', [.5 .5 .5]);
h = plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^2/150000 ,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
l = legend('fast-SLS','gurobi','mosek','$\mathcal{O}(N^{2.5})$','$\mathcal{O}(N^2)$','interpreter','latex');
l.Position = [    0.2327    0.7521    0.1582    0.1322];
xlabel('Horizon length N','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
N_min = [min([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
N_max = [max([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
cpu_min = [min([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
cpu_max = [max([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
axis([N_min, N_max+50, cpu_min/2, cpu_max*2])


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
grid on;
%exportgraphics(gcf,strcat('img/fig1.pdf'),'ContentType','vector');
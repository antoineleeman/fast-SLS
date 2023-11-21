clear all;
close all;
clc;
%addpath(genpath('../casadi_linux'));
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'));
import casadi.*
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

%%

grid_density = 5;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);

timings_N = [];
%%
for nn=5:5:50
    nn
    kkt = KKT_SLS(nn,Q,R,m,Qf); %x0 seems unused %% check if the bo are well reset
    timing = [];
    for ii = 1:length(x1_range)
        for jj = 1:length(x2_range)
            x0 = [x1_range(ii); x2_range(jj)];
            tic
            feasible = kkt.solve(x0);
            time =toc;
            if feasible
                timing = [timing;time];
            end
        end
    end
    timings_N = [timings_N,[nn; mean(timing);std(timing)]];
end
%%

timings_N_yal = [];


for nn=5:10:50
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,m,Qf); %x0 seems unused %% check if the bo are well reset
    timing = [];
    for ii = 1:length(x1_range)
        for jj = 1:length(x2_range)
            x0 = [x1_range(ii); x2_range(jj)];
            tic
            feasible = solver_yalmip.solve(x0);
            time =toc;
            if feasible
                timing = [timing;time];
            end
        end
    end
    timings_N_yal = [timings_N_yal,[nn; mean(timing);std(timing)]];
end

%%
load('long_horizon');
clf;
errorbar(timings_N(1,:), timings_N(2,:), timings_N(3,:),'LineWidth',2);
hold on;
errorbar(timings_N_yal(1,:), timings_N_yal(2,:), timings_N_yal(3,:),'LineWidth',2);
plot(timings_N(1,:), timings_N(1,:).^2/1000 ,'LineWidth',2);
plot(timings_N(1,:), timings_N(1,:).^4/10000 ,'LineWidth',2);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('iSLS','gurobi','$\mathcal{O}(N^2)$','$\mathcal{O}(N^4)$','interpreter','latex');
xlabel('N');
ylabel('computation times');

set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 15]);

exportgraphics(gcf,strcat('fig1.pdf'),'ContentType','vector');

%%
clear all;
close all;
clc;

N = 15;


n_sample = 3;

for mm = 3:3:15
    msd = ChainOfMassSpringDampers(3);
    Q = eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf); 
    kkt = KKT_SLS(N,Q,R,msd,Qf); %x0 seems unused %% check if the bo are well reset

    timing_mm_yalmip = [];
    timing_mm_kkt = [];
    for ii=1:n_sample
        tic
        feasible = solver_yalmip.solve(rand(msd.nx,1));
        time =toc;
        if feasible
            timing_mm_yalmip = [timing_mm_yalmip;time];
        end
    
        tic
        feasible = kkt.solve(rand(msd.nx,1));
        time =toc;
        if feasible
            timings_M_kkt = [timing_mm_kkt;time];
        end
    end

   timings_M_yal = [timings_M_yal,[mm; mean(timing_mm_yalmip);std(timing_mm_yalmip)]];
   timings_M_kkt = [timings_M_kkt,[mm; mean(timings_M_kkt);std(timings_M_kkt)]];


end

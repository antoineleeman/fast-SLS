clear all;
close all;
clc;
%addpath(genpath('../casadi_linux'));
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'));
import casadi.*
m = Integrator();
%x0 = [-2;5];
x0 = [-2.5;5];
N =10;
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

%%

grid_density = 5;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);

timing_explicit = [];

timings_N = [];
%%
for nn=5:5:50
    nn
    kkt = KKT_SLS(nn,Q,R,m,x0,Qf); %x0 seems unused %% check if the bo are well reset
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

for nn=5:5:50
    nn
    solver_sedumi = YALMIP_SLS(nn,Q,R,m,x0,Qf); %x0 seems unused %% check if the bo are well reset
    timing = [];
    for ii = 1:length(x1_range)
        for jj = 1:length(x2_range)
            x0 = [x1_range(ii); x2_range(jj)];
            tic
            feasible = solver_sedumi.solve(x0);
            time =toc;
            if feasible
                timing = [timing;time];
            end
        end
    end
    timings_N = [timings_N,[nn; mean(timing);std(timing)]];
end

%%
clf;
errorbar(timings_N(1,:), timings_N(2,:), timings_N(3,:),'LineWidth',2);
hold on;
plot(timings_N(1,:), timings_N(1,:).^2/1000 ,'LineWidth',2);
plot(timings_N(1,:), timings_N(1,:).^4/10000 ,'LineWidth',2);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('iSLS', '$\mathcal{O}(N^2)$','$\mathcal{O}(N^4)$','interpreter','latex');
xlabel('N');
ylabel('computation times');

%%
x0 = [-10;5];

solver_sedumi = YALMIP_SLS(5,Q,R,m,x0,Qf);
solver_sedumi.solve(x0)



clear all;
close all;
clc;
%%
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

grid_density = 15;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
N = 3;
%
profile on
kkt = KKT_SLS(N,Q,R,m,Qf); %% check if the bo are well reset
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end
%%
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

grid_density = 15;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
%
profile on
kkt = KKT_SLS(15,Q,R,m,Qf); %% check if the bo are well reset
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end


%%
L = 3;
msd = ChainOfMassSpringDampers(L);
Q = eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
grid_density = 5;
% x1_range = linspace(-5,5,grid_density);
% x2_range = linspace(-5,5,grid_density);

timings_N_exact_kkt = [];
%for nn=5:15:100
n_sample = 15;
for nn=5:10:80
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf); %x0 seems unused %% check if the bo are well reset
    timing_kkt = [];

    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
        feasible = kkt.solve(x0);
        time =toc;
        if feasible
            timing_kkt = [timing_kkt;time];
        end
    end
    %     end
    % end
    timings_N_exact_kkt = [timings_N_exact_kkt,[nn; mean(timing_kkt);std(timing_kkt)]];
    
end
save(getUniqueName('timings_N_exact_kkt'),'timings_N_exact_kkt','msd')

%%

timings_N_yal = [];
for nn=5:10:80
    nn
    solver_yalmip = YALMIP_SLS(nn,Q,R,msd,Qf); %x0 seems unused %% check if the bo are well reset
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
    timings_N_yal = [timings_N_yal,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_yal'),'timings_N_yal','msd')

%%
clear all;
close all;


clf;

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];


errorbar(timings_N_exact_kkt(1,:), timings_N_exact_kkt(2,:), timings_N_exact_kkt(3,:),'LineWidth',2,'Color', colors(1,:));
hold on;
errorbar(timings_N_yal(1,:), timings_N_yal(2,:), timings_N_yal(3,:),'LineWidth',2,'Color', colors(3,:));
%plot(timings_N_kkt(1,:), timings_N_kkt(1,:).^2/100 ,'LineWidth',2, 'Color', [.5 .5 .5]);
plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(2.5)/1000 ,'LineWidth',2, 'Color', [.5 .5 .5]);
h = plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^2/150000 ,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
l = legend('iSLS','gurobi','$\mathcal{O}(N^{2.5})$','$\mathcal{O}(N^2)$','interpreter','latex');
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

%%
clear all;
close all;
clc;

N = 15;


n_sample = 3;

timings_M_kkt = [];
for mm = 5:20:150
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
        feasible = kkt.solve(x0);
        time =toc;
        if feasible
            timings_mm_kkt = [timing_mm_kkt;time];
        end
    end

    timings_M_kkt = [timings_M_kkt,[mm; mean(timings_mm_kkt);std(timings_mm_kkt)]];
end

%%
timings_M_yal = [];
N = 15;
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

save('yalm_scalability_nx');
%%
figure(2);
clf;
load('yalm_scalability_nx.mat');
load('kkt_scalability_nx.mat');

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
legend('iSLS','Gurobi','$\mathcal{O}(n_x^3)$','interpreter','latex');
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
exportgraphics(gcf,strcat('img/fig2.pdf'),'ContentType','vector');


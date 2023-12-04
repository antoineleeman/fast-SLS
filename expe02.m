
%%
clear all;
close all;
clc;

%N = 15;
%n_sample = 1;

%%
expe01_init
timings_M_kkt = [];
IT_M_kkt = [];
for mm = 3:10:80
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
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

histogram(IT_M_kkt)
save(getUniqueName('timings_M_kkt'),'timings_M_kkt','IT_M_kkt','msd','N','n_sample')
save('timings_M_kkt.mat','timings_M_kkt','msd','N','n_sample');


%%
expe01_init
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
        x0 =rand(msd.nx,1);
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
expe01_init
timings_M_mosek = [];
for mm = 3:1:6
    mm
    msd = ChainOfMassSpringDampers_actuated(mm);
    Q = 3*eye(msd.nx);
    R = eye(msd.nu);
    Qf = Q;
    solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'mosek'); 
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
N
load('timings_M_gurobi');
N
load('timings_M_kkt');
N


%load('yalm_scalability_nx.mat');
%load('kkt_scalability_nx.mat');
%load('24-Nov-2023_16_30_08__timings_M_yal.mat');
%load('24-Nov-2023_16_32_27__timings_M_kkt.mat');

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];


plot(timings_M_kkt(1,:), timings_M_kkt(2,:),'LineWidth',2,'Color', colors(1,:),'Marker','x');
hold on;
plot(timings_M_gurobi(1,:), timings_M_gurobi(2,:),'LineWidth',2,'Color', colors(3,:),'Marker','+');
plot(timings_M_mosek(1,:), timings_M_mosek(2,:),'LineWidth',2,'Color', colors(4,:),'Marker','+');

% First plot with DisplayName for legend
plot(timings_M_kkt(1,:), timings_M_kkt(1,:).^3/500000, 'LineWidth', 2, 'Linestyle',':', 'DisplayName', '$\mathcal{O}(n_x^3)$', 'Color', [.5 .5 .5]);

% Second plot without legend entry
h = plot(timings_M_mosek(1,:), timings_M_mosek(1,:).^3/50, 'LineWidth', 2, 'Linestyle',':', 'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

% Second plot without legend entry
h = plot(timings_M_gurobi(1,:), timings_M_gurobi(1,:).^3/10, 'LineWidth',2,'Linestyle',':',  'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
legend('fast-SLS','gurobi','mosek','$\mathcal{O}(n_x^3)$','interpreter','latex');
xlabel('Number of states $n_x$','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
% nx_min = [min([timings_M_gurobi(1,:),timings_M_kkt(1,:)])];
% nx_max = [max([timings_M_gurobi(1,:),timings_M_kkt(1,:)])];
% cpu_min = [min([timings_M_gurobi(2,:),timings_M_kkt(2,:)])];
% cpu_max = [max([timings_M_gurobi(2,:),timings_M_kkt(2,:)])];
% axis([nx_min, nx_max, cpu_min/2, cpu_max*2])
axis([timings_M_kkt(1,1), timings_M_kkt(1,end)+50, 0.005, 150])


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
exportgraphics(gcf,strcat('img/fig2.pdf'),'ContentType','vector');
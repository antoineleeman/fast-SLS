
%%
expe01_init
timings_N_exact_kkt = [];

for nn=3:5:300
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf);
    timing_kkt = [];

    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
        [feasible] = kkt.solve(x0);
        time =toc;
        if feasible
            [timing_kkt] = [timing_kkt;time];
        end
    end
    if time> 60
        break
    end
    %     end
    % end
    timings_N_exact_kkt = [timings_N_exact_kkt,[nn; mean(timing_kkt);std(timing_kkt)]];
    
end

save(getUniqueName('timings_N_exact_kkt'),'timings_N_exact_kkt','msd')
%%
expe01_init
timings_N_rti_kkt = [];

for nn=3:5:150
    nn
    kkt = RTI_fast_SLS(nn,Q,R,msd,Qf);
    timing_kkt_rti = [];

    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        tic
            [feasible] = kkt.solve(x0);
        time =toc;
        [timing_kkt_rti] = [timing_kkt_rti;time];
    end
    timings_N_rti_kkt = [timings_N_rti_kkt,[nn; mean(timing_kkt_rti);std(timing_kkt_rti)]];
end

save(getUniqueName('timings_N_rti_kkt'),'timings_N_rti_kkt','msd')


%%

timings_N_gurobi = [];
for nn=3:3:25
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
        if time> 60
        break
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
    if time> 60
        break
    end
    timings_N_mosek = [timings_N_mosek,[nn; mean(timing_yal);std(timing_yal)]];
end
save(getUniqueName('timings_N_mosek'),'timings_N_mosek','msd')
%%
expe01_init
timings_N_nominal = [];

for nn=3:5:300
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf);
    timing_nom = [];

    for ii =1:n_sample
        x0 =rand(msd.nx,1);
        [~, time] = kkt.forward_solve(x0);
        [timing_nom] = [timing_nom;time];
    end
    %     end
    % end
    timings_N_nominal = [timings_N_nominal,[nn; mean(timing_nom);std(timing_nom)]];
    
end
save(getUniqueName('timings_N_nominal'),'timings_N_nominal','msd');



%%
clear all;
close all;
clf;
load('25-Nov-2023_15_38_00__timings_N_gurobi.mat')
msd.nx
load('25-Nov-2023_16_47_49__timings_N_exact_kkt.mat')
msd.nx
load('25-Nov-2023_15_13_02__timings_N_mosek.mat')
msd.nx
%load('25-Nov-2023_16_28_54__timings_N_rti_kkt.mat')
msd.nx
load('25-Nov-2023_16_50_39__timings_N_nominal.mat')
msd.nx


clf;

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
figure(1);
plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(2,:),'LineWidth',2,'Color', colors(1,:));
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(timings_N_gurobi(1,:), timings_N_gurobi(2,:),'LineWidth',2,'Color', colors(2,:));
plot(timings_N_mosek(1,:), timings_N_mosek(2,:),'LineWidth',2,'Color', colors(3,:));
plot(timings_N_nominal(1,:), timings_N_nominal(2,:),'LineWidth',2,'Color', colors(4,:));
%plot(timings_N_rti_kkt(1,:), timings_N_rti_kkt(2,:),'LineWidth',2,'Color', colors(5,:));

plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(2.5)/1000 ,'LineWidth',2, 'Linestyle',':', 'Color', [.5 .5 .5]);
%h = plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^2/150000 ,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');


l = legend('fast-SLS','gurobi','mosek','nominal','$\mathcal{O}(N^{2.5})$','interpreter','latex');
l.Position = [    0.2327    0.7521    0.1582    0.1322];
xlabel('Horizon length N','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
% N_min = [min([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% N_max = [max([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% cpu_min = [min([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
% cpu_max = [max([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
% axis([N_min, N_max+50, cpu_min/2, cpu_max*2])


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
grid on;
%exportgraphics(gcf,strcat('img/fig1.pdf'),'ContentType','vector');
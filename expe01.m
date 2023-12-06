
%%
clear all;
close all;
clc;


%%
expe01_init
timings_N_exact_kkt = [];
IT = [];
for nn=3:10:120
    nn
    kkt = KKT_SLS(nn,Q,R,msd,Qf);
    timing_kkt = [];

    for ii =1:n_sample
        x0 =4*rand(msd.nx,1)-2;
        tic
        [feasible, it ] = kkt.solve(x0);
        time =toc;
        if feasible
            [timing_kkt] = [timing_kkt;time];
            IT = [IT,it];
        end
    end

    timings_N_exact_kkt = [timings_N_exact_kkt,[nn; mean(timing_kkt);std(timing_kkt)]];
    
end
histogram(IT)


%save(getUniqueName('timings_N_exact_kkt'),'timings_N_exact_kkt','msd')
%save('fast-sls-N.mat','timings_N_exact_kkt','msd');
%%
expe01_init
timings_N_rti_kkt = [];

for nn=3:30:120
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
save('rti-fast-sls-N.mat','timings_N_rti_kkt','msd')


%%
expe01_init;
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
expe01_init
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

%%
expe01_init
timings_N_nominal = [];

for nn=3:15:150
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
save('nominal-N.mat','timings_N_nominal','msd');



%%
clear all;
close all;
clf;
load('gurobi-sls-N.mat')
msd.nx
load('fast-sls-N.mat') % use this: 05-Dec-2023_14_17_22__timings_N_exact_kkt
msd.nx
load('mosek-sls-N.mat')
msd.nx
%load('rti-fast-sls-N.mat')
msd.nx
%load('nominal-N.mat')
msd.nx


clf;

colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
figure(1);
errorbar(timings_N_exact_kkt(1,:), timings_N_exact_kkt(2,:),timings_N_exact_kkt(3,:),'LineWidth',2,'Color', colors(1,:),'Marker','x');
hold on;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(timings_N_gurobi(1,:), timings_N_gurobi(2,:),'LineWidth',2,'Color', colors(3,:),'Marker','+');
plot(timings_N_mosek(1,:), timings_N_mosek(2,:),'LineWidth',2,'Color', colors(4,:),'Marker','+');
%plot(timings_N_nominal(1,:), timings_N_nominal(2,:),'LineWidth',2,'Color', colors(4,:));
%plot(timings_N_rti_kkt(1,:), timings_N_rti_kkt(2,:),'LineWidth',2,'Color', colors(3,:));

plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(1)/5500 ,'LineWidth',2, 'Linestyle',':', 'Color', [.5 .5 .5]);

h = plot(timings_N_mosek(1,:), timings_N_mosek(1,:).^3/500 ,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
h = plot(timings_N_gurobi(1,:), timings_N_gurobi(1,:).^2/10 ,'LineWidth',2, 'Linestyle','-.','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');



l = legend('fast-SLS','gurobi','mosek','$\mathcal{O}(N)$','$\mathcal{O}(N^{3})$','$\mathcal{O}(N^{2})$','interpreter','latex');
l.Position = [     0.6506    0.6926    0.1582    0.1322
];
xlabel('Horizon length N','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
% N_min = [min([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% N_max = [max([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% cpu_min = [min([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
% cpu_max = [max([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
%axis([timings_N_exact_kkt(1,1), timings_N_exact_kkt(1,end)+50, 0.005, 50])


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 10]);
exportgraphics(gcf,strcat('img/fig1.pdf'),'ContentType','vector');
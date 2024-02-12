
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
fontsize = 12;
colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
figure(1);
subplot(1,2,1);
errorbar(timings_N_exact_kkt(1,:), timings_N_exact_kkt(2,:),timings_N_exact_kkt(3,:),'LineWidth',2,'Color', colors(1,:));
hold on;
errorbar(timings_N_exact_ff(1,:), timings_N_exact_ff(2,:), timings_N_exact_ff(3,:),'LineStyle','--','LineWidth',2,'Color', colors(1,:));

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
plot(timings_N_gurobi(1,:), timings_N_gurobi(2,:),'LineWidth',2,'Color', colors(3,:),'Marker','+');
plot(timings_N_mosek(1,:), timings_N_mosek(2,:),'LineWidth',2,'Color', colors(4,:),'Marker','+');
%plot(timings_N_nominal(1,:), timings_N_nominal(2,:),'LineWidth',2,'Color', colors(4,:));
%plot(timings_N_rti_kkt(1,:), timings_N_rti_kkt(2,:),'LineWidth',2,'Color', colors(3,:));

plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(2)/55000 ,'LineWidth',2, 'Linestyle','-.', 'Color', [.5 .5 .5]);

h = plot(timings_N_mosek(1,:), timings_N_mosek(1,:).^3/100 /2,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
h = plot(timings_N_gurobi(1,:), timings_N_gurobi(1,:).^2/2.5/2 ,'LineWidth',2, 'Linestyle','-.','Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');



l1 = legend('fast-SLS','fast-SLS: riccati','gurobi','mosek','$\mathcal{O}(N^2)$','$\mathcal{O}(N^{3})$','interpreter','latex','FontSize',fontsize);
l1.Position = [0.3092 0.5775 0.1621 0.3283];
set(l1, 'Box', 'off', 'Color', 'none');
set(gca,'FontSize',12);

xlabel('Horizon length N','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
% N_min = [min([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% N_max = [max([timings_N_yal(1,:),timings_N_exact_kkt(1,:)])];
% cpu_min = [min([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
% cpu_max = [max([timings_N_yal(2,:),timings_N_exact_kkt(2,:)])];
axis([timings_N_exact_kkt(1,1), timings_N_exact_kkt(1,end)+10, 0.001, 100])


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
subplot(1,2,2);

errorbar(timings_M_kkt(1,:), timings_M_kkt(2,:), timings_M_kkt(3,:),'LineWidth',2,'Color', colors(1,:));

hold on;

errorbar(timings_M_kkt_ff(1,:), timings_M_kkt_ff(2,:), timings_M_kkt_ff(3,:),'Linestyle','--','LineWidth',2,'Color', colors(1,:));


plot(timings_M_gurobi(1,:), timings_M_gurobi(2,:),'LineWidth',2,'Color', colors(3,:),'Marker','+');
plot(timings_M_mosek(1,:), timings_M_mosek(2,:),'LineWidth',2,'Color', colors(4,:),'Marker','+');

% First plot with DisplayName for legend
plot(timings_M_kkt(1,:), timings_M_kkt(1,:).^3/300000, 'LineWidth', 2, 'Linestyle',':', 'DisplayName', '$\mathcal{O}(n_x^3)$', 'Color', [.5 .5 .5]);

% Second plot without legend entry
%h = plot(timings_M_mosek(1,:), timings_M_mosek(1,:).^3/50, 'LineWidth', 2, 'Linestyle',':', 'Color', [.5 .5 .5]);
%set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

% Second plot without legend entry
h = plot(timings_M_gurobi(1,:), timings_M_gurobi(1,:).^3/10, 'LineWidth',2,'Linestyle',':',  'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
l2 = legend('fast-SLS','fast-SLS: riccati','gurobi','mosek','$\mathcal{O}(n_x^3)$','interpreter','latex','FontSize',fontsize);
set(l2, 'Box', 'off', 'Color', 'none');

xlabel('Number of states $n_x$','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
% nx_min = [min([timings_M_gurobi(1,:),timings_M_kkt(1,:)])];
% nx_max = [max([timings_M_gurobi(1,:),timings_M_kkt(1,:)])];
% cpu_min = [min([timings_M_gurobi(2,:),timings_M_kkt(2,:)])];
% cpu_max = [max([timings_M_gurobi(2,:),timings_M_kkt(2,:)])];
% axis([nx_min, nx_max, cpu_min/2, cpu_max*2])
axis([timings_M_kkt(1,1), timings_M_kkt(1,end)+10, 0.005, 150])

set(gca,'FontSize',12);
set(gcf,'units','centimeters','Position', [0 0 16.1*2 10]);
exportgraphics(gcf,strcat('img/fig4.pdf'),'ContentType','vector');

%%
%%
clear all
close all
clf
%load('26-Nov-2023_17_23_21__it_kkt.mat')
%load('28-Nov-2023_13_23_07__it_kkt.mat')

load('06-Dec-2023_17_34_03__it_kkt.mat');
colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
h = histogram(it_kkt, 'FaceColor', colors(1,:));
xlabel('Number of iterations','interpreter','latex');

% Total number of simulations
numSimulations = numel(it_kkt);

% Specific percentages you want to display
specificPercentages = [0.1, 1, 10, 100];

% Convert these percentages to corresponding counts
specificCounts = specificPercentages * numSimulations / 100;

% Set logarithmic scale for y-axis
set(gca, 'YScale', 'log');

% Set the y-axis to these specific counts and label them with the corresponding percentages
set(gca, 'YTick', specificCounts, 'YTickLabel', strcat(num2str(specificPercentages'), '%'));

ylabel('Percentage of simulations','interpreter','latex');
grid on;


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 6]);
exportgraphics(gcf,strcat('img/fig3.pdf'),'ContentType','vector');


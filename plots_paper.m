% File: plots_paper.m
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
clf;
load('data/gurobi-sls-N.mat')
load('data/fast-sls-N.mat')
load('data/mosek-sls-N.mat')

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
plot(timings_N_exact_kkt(1,:), timings_N_exact_kkt(1,:).^(2)/55000 ,'LineWidth',2, 'Linestyle','-.', 'Color', [.5 .5 .5]);

h = plot(timings_N_mosek(1,:), timings_N_mosek(1,:).^3/100 /2,'LineWidth',2, 'Linestyle','--','Color', [.5 .5 .5]);
h = plot(timings_N_gurobi(1,:), timings_N_gurobi(1,:).^2/2.5/2 ,'LineWidth',2, 'Linestyle','-.','Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');



l1 = legend('fast-SLS','fast-SLS: riccati','gurobi','mosek','$\mathcal{O}(N^2)$','$\mathcal{O}(N^{3})$','interpreter','latex','FontSize',fontsize);
l1.Position = [0.3092 0.5775 0.1621 0.3283];
set(l1, 'Box', 'off', 'Color', 'none');
set(gca,'FontSize',12);

xlabel('Horizon length N','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;
axis([timings_N_exact_kkt(1,1), timings_N_exact_kkt(1,end)+10, 0.001, 100])

load('data/timings_M_mosek');
load('data/timings_M_gurobi');
load('data/timings_M_kkt');

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
plot(timings_M_kkt(1,:), timings_M_kkt(1,:).^3/300000, 'LineWidth', 2, 'Linestyle',':', 'DisplayName', '$\mathcal{O}(n_x^3)$', 'Color', [.5 .5 .5]);
h = plot(timings_M_gurobi(1,:), timings_M_gurobi(1,:).^3/10, 'LineWidth',2,'Linestyle',':',  'Color', [.5 .5 .5]);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');


set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
l2 = legend('fast-SLS','fast-SLS: riccati','gurobi','mosek','$\mathcal{O}(n_x^3)$','interpreter','latex','FontSize',fontsize);
set(l2, 'Box', 'off', 'Color', 'none');

xlabel('Number of states $n_x$','interpreter','latex');
ylabel('Computation time [s]','interpreter','latex');
grid on;

axis([timings_M_kkt(1,1), timings_M_kkt(1,end)+10, 0.005, 150])
xticks([10^1 70]);
set(gca,'FontSize',12);
set(gcf,'units','centimeters','Position', [0 0 16.1*2 10]);
exportgraphics(gcf,strcat('img/fig4.pdf'),'ContentType','vector');
%%
clear all
close all
clf

load('data/06-Dec-2023_17_34_03__it_kkt.mat');
colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
h = histogram(it_kkt, 'FaceColor', colors(1,:));
xlabel('Number of iterations','interpreter','latex');

numSimulations = numel(it_kkt);
specificPercentages = [0.1, 1, 10, 100];
specificCounts = specificPercentages * numSimulations / 100;
set(gca, 'YScale', 'log');
set(gca, 'YTick', specificCounts, 'YTickLabel', strcat(num2str(specificPercentages'), '%'));
ylabel('Percentage of simulations','interpreter','latex');
grid on;


set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 6]);
exportgraphics(gcf,strcat('img/fig3.pdf'),'ContentType','vector');


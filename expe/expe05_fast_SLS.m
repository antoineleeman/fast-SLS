% File: expe05_fast_SLS.m
% Author: Antoine Leeman (aleeman(at)ethz(dot)ch)
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
N = 20;
m = Integrator();
x_max = 5;
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = 100*Q;

kkt = KKT_SLS(N,Q,R,m,Qf);
x0 = [-3; 5];
[feasible,it, ~, ~, ~, ~,it_data] = kkt.solve(x0);
        
figure(1);
clf;
colors =[    0.2670    0.0049    0.3294
    0.2673    0.2259    0.5132
    0.1906    0.4071    0.5561
    0.1282    0.5651    0.5509
    0.2080    0.7187    0.4729
    0.5604    0.8413    0.2661
    0.9932    0.9062    0.1439
    ];

color_MPSF = colors(1,:);
color_Se = colors(2,:);
color_SL_MPSF = colors(3,:);
color_max_RPI = colors(4,:);
color_X = colors(5,:);
color_max_RCI = colors(6,:);

fontsize = 18;
hold on;
xlabel('State $x_1$','Interpreter','latex', FontSize=fontsize);
ylabel('State $x_2$','Interpreter','latex', FontSize=fontsize);
grid on;
axis equal;
set(gcf,'units','centimeters','Position', [0 0 25 15]);

rectangle('Position',[-x_max -x_max 2*x_max 2*x_max],'EdgeColor','k','LineStyle','-','Linewidth',2);
n_iter = it_data(1).n_iter;
cm = viridis(n_iter);
XLim =[   -2.2949    6.2901];
YLim = [   -1.5878    5.1832];

% text(5, 0.5, 'Constraint $\mathcal{X}$', 'FontSize',fontsize,'Color', color_X,Interpreter='latex');

for jj = 1 : n_iter-1
    jj
    X = full(it_data(jj).x);
    plot(X(1,:), X(2,:),'Color', cm(jj,:),'LineStyle','-', 'LineWidth', 2);
    str1 = 'traj';
    str1b = 'tube';
    str2 = num2str(jj);
    str3 = '.png';
    xlim(XLim); % Set x-axis limits
    ylim(YLim); % Set y-axis limits
    if jj==2
        bo = it_data(jj).bo_j(:,5);
        rectangle('Position',[-x_max+bo(5) -x_max+bo(3) 2*x_max-2*bo(5) 2*x_max-2*bo(3)],'EdgeColor','k','LineStyle','--','Linewidth',2);
    end
    saveas(gcf, [str1, str2, str3]);
    %for ii = 1:N-1%% hard-coded value!
    for ii = 1:5%% hard-coded value!
        ii
        ellipse_points = MinkowskiSumEllipsoids_bis(it_data(jj).Phi_x(ii,:));
        alpha_value = 0.25;  % Example transparency level (0 = fully transparent, 1 = fully opaque)

        % Plot the patch
        fill(ellipse_points(1, :) + X(1, ii+1), ellipse_points(2, :) + X(2, ii+1), cm(jj,:), 'EdgeColor', 'none', 'FaceAlpha', alpha_value);

        % plot(ellipse_points(1, :) + X(1,ii+1), ellipse_points(2, :)+ X(2,ii+1),'Color', cm(jj,:), 'LineWidth', 2);
        hold on
    end


    xlim(XLim); % Set x-axis limits
    ylim(YLim); % Set y-axis limits

    saveas(gcf, [str1b, str2, str3]);

end

%%
clear all;
close all;
clc;

L = 25;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 1;
N=25;



kkt = KKT_SLS(N,Q,R,msd,Qf);
it_kkt = [];

alpha_delta = zeros(kkt.MAX_ITER,10);
ii=1;
for alpha =0:0.1:1

    ii
    %x0 =4*rand(msd.nx,1)-2;
    x0 = alpha*2.55*ones(msd.nx,1);
    %x0 = -0.99*3*[ones(msd.nx,1)];
    %x0 = kron(ones(msd.nx/2,1),[1;0]);
    %x0= [0;1;zeros(msd.nx-2,1)];
    tic
    [feasible,it,delta] = kkt.solve(x0);
    time =toc;

    alpha_delta(:,ii) = cell2mat(delta);

    ii = ii+1;
end
disp('percentage solved');
length(it_kkt)/n_sample
%
%save(getUniqueName('it_kkt'),'it_kkt','msd');
%
%clear all
%close all
%clf
%load('26-Nov-2023_17_23_21__it_kkt.mat')
%load('28-Nov-2023_13_23_07__it_kkt.mat')
colors = [0.0504    0.0298    0.5280
    0.4934    0.0115    0.6580
    0.7964    0.2780    0.4713
    0.9722    0.5817    0.2541
    0.9400    0.9752    0.1313];
h = histogram(it_kkt, 'FaceColor', colors(1,:));
xlabel('Number of iterations','interpreter','latex');
ylabel('Number of simulations','interpreter','latex');
grid on;

set(gca,'FontSize',10);
set(gcf,'units','centimeters','Position', [0 0 15 6]);
%exportgraphics(gcf,strcat('img/fig4.pdf'),'ContentType','vector');
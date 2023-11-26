
%%
clear all;
close all;
clc;

%            obj.E = 0.05*eye(obj.nw);

%            u_max = 1;
%            x_max = 3;


L = 30;
msd = ChainOfMassSpringDampers(L);
Q = 100*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 100;
N=25;



kkt = KKT_SLS(N,Q,R,msd,Qf);
it_kkt = [];

for ii =1:n_sample
    ii
    x0 =rand(msd.nx,1);
    tic
    [feasible,it] = kkt.solve(x0);
    time =toc;
    if feasible
        [it_kkt] = [it_kkt;it];
    end
end

%save(getUniqueName('it_kkt'),'it_kkt','msd');
%
h = histogram(it_kkt);
xlabel('Number of iterations','interpreter','latex');
ylabel('Number of simulations','interpreter','latex');
grid on;

% set(gca,'FontSize',10);
% set(gcf,'units','centimeters','Position', [0 0 15 6]);
% exportgraphics(gcf,strcat('img/fig3.pdf'),'ContentType','vector');
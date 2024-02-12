
%%
clear all;
close all;
clc;

L = 4;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
N=10;



kkt = KKT_SLS(N,Q,R,msd,Qf);
solver_yalmip = YALMIP_SLS(N,Q,R,msd,Qf,'mosek');

x0 =4*rand(msd.nx,1)-2;
[feasible,it] = kkt.solve(x0);
feasible = solver_yalmip.solve(x0);


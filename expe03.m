
%%
clear all;
close all;
clc;

L = 6;
msd = ChainOfMassSpringDampers(L);
Q = eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 1000;
IT = [];
N=30;

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

histogram(it_kkt)
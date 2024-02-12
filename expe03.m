
%%
clear all;
close all;
clc;

L = 25;
msd = ChainOfMassSpringDampers_actuated(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 1000;
N=25;



kkt = KKT_SLS(N,Q,R,msd,Qf);
it_kkt = [];

for ii =1:n_sample
    ii
    x0 =4*rand(msd.nx,1)-2;
    tic
    [feasible,it] = kkt.solve(x0);
    time =toc;
    if feasible
        [it_kkt] = [it_kkt;it];
    end
end
disp('percentage solved');
length(it_kkt)/n_sample

save(getUniqueName('it_kkt'),'it_kkt','msd');

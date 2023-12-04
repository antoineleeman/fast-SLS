
%%
clear all;
close all;
clc;

%            obj.E = 0.05*eye(obj.nw);

%            u_max = 1;
%            x_max = 3;
L = 30;
msd = ChainOfMassSpringDampers(L);
Q = 3*eye(msd.nx);
R = eye(msd.nu);
Qf = Q;
n_sample = 10;
N=20;
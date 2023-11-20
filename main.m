clear all;
close all;
clc;
%addpath(genpath('../casadi_linux'));
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'));
import casadi.*
m = Integrator();
%x0 = [-2;5];
x0 = [-5;5];
N =10;
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

kkt = KKT_SLS(N,Q,R,m, x0,Qf); %x0 seems unused
kkt.solve();
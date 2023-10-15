v%%
clear all;
close all;
clc;

%addpath('util')
%addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'))

% Import CasADi toolbox
import casadi.*

%try a simple HPIPM

y_bar = MX.sym('y_bar',length(y),1);
Dy = MX.sym('Dy',length(y),1);

f_qp = dfdy_fun(y_bar)* Dy + 0.5*Dy'*H_tilde*Dy;
g_fun = casadi.Function('g_fun',{y},{g});

dGdy_fun = casadi.Function('dGdy_fun',{y},{dGdy});

g_lin = g_fun(y_bar) + dGdy_fun(y_bar) * Dy;

qp = struct('x',Dy, 'f',f_qp, 'g',g_lin,'p',[y_bar]);

%%
addpath(genpath('../casadi_linux'));
import casadi.*

A=[1, 0.2, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0;
 -0.1, 0.4, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
 0.3, 0.2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
 2, 0, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 1, 1, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
 0, 0, 0, 1, 4, 2, 1, 0.3, -1, 0, 0, 0; 
 0, 0, 0, 3, 1, 0, 1, 0.2, 0, -1, 0, 0; 
 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, -1; 
 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 1, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3];

A = DM(sparse(A));

H=[7, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
 0, 7, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
 0.2, 0.3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
 0, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0; 
 0, 0, 0, 0, 2, 0.1, 0, 0.7, 0, 0, 0, 0; 
 0, 0, 0, 0, 0.1, 1, 0, 1, 0, 0, 0, 0; 
 0, 0, 0, 0, 0, 0, 1, 0.1, 0, 0, 0, 0; 
 0, 0, 0, 1, 0.7, 1, 0.1, 2, 0, 0, 0, 0; 
 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 1, 0; 
 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0; 
 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 4, 0; 
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9];

H = DM(sparse(H));

nx = [2,3,2,1];
nu = [1, 2,1];
ng = [2, 1, 1, 1];
N = 3;
bo = MX.sym('bo',1,1);

options = struct;
% options.hpipm.iter_max = 100;
% options.hpipm.res_g_max = 1e-10;
% options.hpipm.res_b_max = 1e-10;
% options.hpipm.res_d_max = 1e-10;
% options.hpipm.res_m_max = 1e-10;

solver = conic('solver', 'qpoases', struct('a',A.sparsity(), 'h', H.sparsity()),options);

g = [1;1;0.2;0.4;1;0.5;0.3;1;0.6;1;1;0.7];
lbg =[0;0;0;-2;-2;0;0;-2;0;-2;-2];
ubg = [0;0;0;2;2;0;0;2;0;2;2];
lbx = [0.5;0.2;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
ubx = [0.5;0.2; 1; 1; 1; 1; 1; 1; 1; 1; 1; bo];

tic
% sol = solver('a',A,'h',H,'lba',lbg,'uba',ubg,'g',g,'lbx',lbx,'ubx',ubx);
%sol = solver('a',A,'h',H,'g',g,'lbx',lbx,'ubx',ubx);
sol = solver('a',A,'h',H,'lba',lbg,'uba',ubg,'g',g);
toc
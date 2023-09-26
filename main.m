%%
clear all;
close all;
clc;

addpath('util')
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'))

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


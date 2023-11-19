clear all;
close all;
clc;
%addpath(genpath('../casadi_linux'));
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'));
import casadi.*
m = Integrator();
x0 = [-5;5];
N =10;
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

kkt = KKT_SLS(N,Q,R,m, x0,Qf); %x0 seems unused

figure(1);
clf;
hold on;
rectangle('Position',[-5 -5 10 10]);
axis equal;

figure(2);
clf;
MAX_ITER = 5;

for ii=1:MAX_ITER

    
    
    [kkt, x_bar, u_bar, lambda_bar, mu_bar] = kkt.forward_solve(x0);
    
    % lambda_bar: dynamics
    % mu_bar: constraints % not sure about their value, are the constraints on
    % x really not active?
    x_bar_plot = full(x_bar);
    u_bar_plot = full(u_bar);

    figure(1);
    plot(x_bar_plot(1,:)+ii,x_bar_plot(2,:),'.-');
    hold on;

    figure(2);
    hold on;
     plot(u_bar_plot);
    
    
    kkt = kkt.update_cost_tube();
    
    [kkt, K] = kkt.backward_solve();

    [kkt,beta] = kkt.udpate_backoff();

end

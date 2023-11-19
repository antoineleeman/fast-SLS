clear all;
close all;
clc;
%addpath(genpath('../casadi_linux'));
addpath(genpath('../casadi-3.6.3-osx64-matlab2018b'));
import casadi.*
m = Integrator();
x0 = [-5;5];

kkt = KKT_SLS(5,eye(m.nx),eye(m.nu),m, x0,eye(m.nx)); %x0 seems unused

figure(1);
hold on;
for ii=1:5
    ii
    [kkt, x_bar, u_bar, lambda_bar, mu_bar] = kkt.forward_solve(x0);
    
    % lambda_bar: dynamics
    % mu_bar: constraints % not sure about their value, are the constraints on
    % x really not active?
    x_bar_plot = full(x_bar);
    u_bar_plot = full(u_bar);

    plot(x_bar_plot(1,:),x_bar_plot(2,:),'.-');
    % figure(2);
    % hold on;
    % plot(u_bar_plot);
    
    
    kkt.update_cost_tube();
    
    [kkt, K] = kkt.backward_solve();
    [kkt,beta] = kkt.udpate_backoff();

end

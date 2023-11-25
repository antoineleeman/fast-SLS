clear all;
close all;
clc;
%%
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

grid_density = 15;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
N = 50;
%
profile on
kkt = KKT_SLS(N,Q,R,m,Qf); %% check if the bo are well reset
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end
%%
m = Integrator();
Q = eye(m.nx);
R = 100*eye(m.nu);
Qf = Q;

grid_density = 15;
x1_range = linspace(-5,5,grid_density);
x2_range = linspace(-5,5,grid_density);
timings_N = [];
%
profile on
kkt = KKT_SLS(15,Q,R,m,Qf); %% check if the bo are well reset
IT = [];
for ii = 1:length(x1_range)
    for jj = 1:length(x2_range)
        x0 = [x1_range(ii); x2_range(jj)];
        [feasible,it] = kkt.solve(x0);
        if feasible
            IT = [IT;it];
        end
    end
end






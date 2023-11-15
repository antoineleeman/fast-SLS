% -----------------------------------------------------------------------------
% File: Integrator.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 15th May 2023
% License: MIT
% Reference:
%{
@inproceedings{leeman2023a,
  title = {Predictive safety filter using system level synthesis},
  year = {2023},
  booktitle = {Proceedings of the 5th Annual Learning for Dynamics and Control Conference, PMLR},
  volume = {211},
  author={Leeman, Antoine P. and K{\"o}hler, Johannes and Benanni, Samir and Zeilinger, Melanie N.},
  pages = {1180-1192},
  doi = {10.3929/ethz-b-000615512}
}
%}
% Link: https://arxiv.org/abs/2212.02111
% -----------------------------------------------------------------------------
classdef Integrator

    properties
        A; %dyn, ct or dt??
        B; %input
        C; % constraint x
        d;
        Cf; % terminal constraint
        df;
        E; % noise matrix
        nx; % number of state variables
        nu; % number of input variables
        nw; % size of the disturbance
        ni; % size of the constraints
        ni_x; %number of terminal constraint
        dt; % time step
    end

    methods
        function obj = Integrator()
            obj.nx=2;
            obj.nu=1;
            obj.ni = 6;
            obj.ni_x =4;
            obj.nw = 2;
            obj.dt =1;
            obj.E = 0.01*eye(obj.nx);
            obj.A = [1,1;...
                0,1];
            obj.B = [0.5;1];
            obj.C = [[1,0;...
                -1,0;...
                0,1;...
                0,-1;
                zeros(2)],[zeros(4,1);
                1;
                -1]];
            x_max = 5;
            u_max = 3;

            obj.d = [x_max;x_max;x_max;x_max;u_max ;u_max];
            obj.Cf = [1,0;...
                -1,0;...
                0,1;...
                0,-1;];
            obj.df = [x_max;x_max;x_max;x_max];
        end

    end
end


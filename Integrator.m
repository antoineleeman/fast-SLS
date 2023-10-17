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
classdef Integrator < DynamicalSystem

    properties
                
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
        end
        function dt = ode(obj,x,u) % equation of motion of the dynamical system in continuous time (x_dot = ode(x,u) )
            A = [1,1;...
                0,1];
            B = [0.5;1];
            dt = A*x+B*u;
        end

        function x_p = ddyn(obj,x,u) %discretization of the dynamical system
            x_p = obj.ode(x,u); %ode already in discrete time
        end

        function [g,f] = cons(obj,x,u)
        % add an if statement for linear constraints
            x_max = 5;
            u_max = 3;
            C = [1,0;...
                -1,0;...
                0,1;...
                0,-1;
                zeros(2)];
            D = [zeros(4,1);
                1;
                -1];
            f =[x_max;x_max;x_max;x_max;u_max ;u_max];
            g = C*x+D*u-f;
        end
        function [g,f] = cons_f(obj,x,u)
            x_max = 5;
            C = [1,0;...
                -1,0;...
                0,1;...
                0,-1];
            f =[x_max;x_max;x_max;x_max];
            g = C*x-f;
        end


    end
end


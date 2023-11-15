% -----------------------------------------------------------------------------
% File: LinearSystem.m
% Author: Antoine Leeman
% Date: 15th November 2023
% License: MIT
% Description: An abstract class that defines a generic linear dynamical system.
% -----------------------------------------------------------------------------

classdef (Abstract) LinearSystem < DynamicalSystem
    properties
        Ass; %dyn, ct
        Bss; %input
        Css; % constraint x
        dss;
        Cfss; % terminal constraint
        dfss;
        Ess; % noise matrix
    end

    methods
        % Constructor
        function obj = LinearSystem(A,B,C,D,E,dt)
            if nargin > 0
                obj.Ass = A;
                obj.Bss = B;
                obj.Css = C;
                obj.Ess = E;
                [obj.nx, obj.nu] = size(B);
                [~,obj.nw] = size(E);
                [obj.ni,~] = size(C);
                [obj.ni_x,~] = size(Cf);
            end
        end
    end

    methods
        function dt = ode(obj,x,u) % equation of motion of the dynamical system in continuous time (x_dot = ode(x,u) )
            dt = obj.A*x+obj.B*u;
        end

        function A = A(obj,x,u)
            A = obj.A;
        end

    end

end

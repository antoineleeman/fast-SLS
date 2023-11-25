% -----------------------------------------------------------------------------
% File: LinearSystem.m
% Author: Antoine Leeman
% Date: 15th November 2023
% License: MIT
% Description: An abstract class that defines a generic linear dynamical system.
% -----------------------------------------------------------------------------

classdef (Abstract) LinearSystem < DynamicalSystem
    properties
        A; %dyn, ct
        B; %input
        C; % constraint x
        d;
        Cf; % terminal constraint
        df;
        E; % noise matrix
    end

    methods
        % Constructor
        function obj = LinearSystem(A,B,C,E,dt)
            if nargin > 0
                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.E = E;
                obj.dt =dt;
                [obj.nx, obj.nu] = size(B);
                [~,obj.nw] = size(E);
                [obj.ni,~] = size(C);
                [obj.ni_x,~] = size(Cf);
            end
        end

        function x_p = ddyn(obj,x,u)
            x_p = obj.A*x+obj.B*u;
        end
    end


end

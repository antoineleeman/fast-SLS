% -----------------------------------------------------------------------------
% File: DynamicalSystem.m
% Author: Antoine Leeman
% Date: 12th October 2023
% License: MIT
% Description: An abstract class that defines a generic dynamical system.
% -----------------------------------------------------------------------------

classdef (Abstract) DynamicalSystem
    properties
        nx; % number of state variables
        nu; % number of input variables
        nw; % size of the disturbance
        ni; % size of the constraints
        x0; % initial state
        dt; % time step
        parameters;
        integrator;
        linear;
    end

    methods (Abstract)
        ode(obj,x,u,w);%continuous time
        setLinear();
    end

    methods
        % Constructor
        function obj = DynamicalSystem(nx, nu, dt)
            if nargin > 0
                obj.nx = nx;
                obj.nu = nu;
                obj.nw = nw;
                obj.dt = dt;
                obj.linear = obj.setLinear();
            end
        end
        function x_p = ddyn(obj,x,u) %discretization of the dynamical system

            integrator = obj.integrator;
            h = obj.dt;
            switch integrator
                case 'single'
                    x_p = x + h*ode(obj,x,u);
                case 'multi'
                    step = 10;
                    for i = 1:step
                        x = x + h/step*ode(obj,x,u);
                    end
                    x_p = x;
                case 'rk4'
                    k_1 = ode(obj,x,u);
                    k_2 = ode(obj,x+0.5*h*k_1,u);
                    k_3 = ode(obj,x+0.5*h*k_2,u);
                    k_4 = ode(obj,x+h*k_3,u);
                    x_p = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
                otherwise
                    error('unrecognised integrator');
            end
        end

        function A = A(obj,x,u) %state matric of the discrete time linearized dynamics
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            A = casadi.jacobian(obj.ddyn(x_fun,u_fun), x_fun);
            A_fun = casadi.Function('A_fun',{var_fun},{A});
            A = A_fun([x;u]);
        end

        function B = B(obj,x,u) %input matrix of the discrete time linearized dynamics
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            B = casadi.jacobian(obj.ddyn(x_fun,u_fun), u_fun);
            B_fun = casadi.Function('B_fun',{var_fun},{B});
            B = B_fun([x;u]);
        end

        function G = G(obj,x,u)
            %TODO
        end
        % add constraints

    end

end

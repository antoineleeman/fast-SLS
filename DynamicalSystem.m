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
        cons(obj,x,u);
    end

    methods
        % Constructor
        function obj = DynamicalSystem(nx, nu,nw, dt)
            if nargin > 0
                obj.nx = nx;
                obj.nu = nu;
                obj.nw = nw;
                obj.dt = dt;
                parameters = struct;
            end
        end
        function x_p = ddyn(obj,x,u) %discretization of the dynamical system

            integrator = 'rk4';
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
            import casadi.*
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            A = jacobian(obj.ddyn(x_fun,u_fun), x_fun);
            A_fun = casadi.Function('A_fun',{var_fun},{A});
            A = A_fun([x;u]);
        end
        function A_sparsity = A_S(obj)
            x = casadi.SX.sym('x',obj.nx);
            u = casadi.SX.sym('u',obj.nu);
            A_sparsity = obj.A(x_fun,u_fun).sparsity();
        end

        function B = B(obj,x,u) %input matrix of the discrete time linearized dynamics
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            B = casadi.jacobian(obj.ddyn(x_fun,u_fun), u_fun);
            B_fun = casadi.Function('B_fun',{var_fun},{B});
            B = B_fun([x;u]);
        end
        function B_sparsity = B_S(obj)
            x = casadi.SX.sym('x',obj.nx);
            u = casadi.SX.sym('u',obj.nu);
            B_sparsity = obj.B(x_fun,u_fun).sparsity();
        end

        function C = C(obj,x,u)
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            C = casadi.jacobian(obj.cons(x_fun,u_fun), x_fun);
            C_fun = casadi.Function('C_fun',{var_fun},{C});
            C = C_fun([x;u]);
        end
        function C_sparsity = C_S(obj)
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            C_sparsity = obj.C(x_fun,u_fun).sparsity();
        end

        function D = D(obj,x,u)
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            D = casadi.jacobian(obj.cons(x_fun,u_fun), U_fun);
            D_fun = casadi.Function('C_fun',{var_fun},{D});
            D = D_fun([x;u]);
        end
        
        function D_sparsity = D_S(obj)
            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            D_sparsity = obj.C(x_fun,u_fun).sparsity();
        end
        function G = G(obj,x,u) % derivative wrt w
            %TODO
        end
        % add constraints

    end

end

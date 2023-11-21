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
        ni_x; %number of terminal constraint
        dt; % time step
        parameters;
        integrator;
        linear; %boolean
    end

    % methods (Abstract)
    %     ode(obj,x,u,w); %continuous time
    %     cons(obj,x,u); %constraint
    %     cons_f(obj,x);
    % end

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


        function B = B(obj,x,u) %input matrix of the discrete time linearized dynamics
                        import casadi.*

            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            B = jacobian(obj.ddyn(x_fun,u_fun), u_fun);
            B_fun = casadi.Function('B_fun',{var_fun},{B});
            B = B_fun([x;u]);
        end



        function C = C(obj,x,u)
            import casadi.*

            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            C = jacobian(obj.cons(x_fun,u_fun), x_fun);
            C_fun = casadi.Function('C_fun',{var_fun},{C});
            C = C_fun([x;u]);
        end

        function Cf = Cf(obj,x) %final constraint
            import casadi.*

            x_fun = casadi.SX.sym('x',obj.nx);
            Cf = jacobian(obj.cons_f(x_fun), x_fun);
            C_fun = casadi.Function('C_fun',{x_fun},{Cf});
            Cf = C_fun(x);
        end


        function D = D(obj,x,u)
            import casadi.*

            x_fun = casadi.SX.sym('x',obj.nx);
            u_fun = casadi.SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            D = jacobian(obj.cons(x_fun,u_fun), u_fun);
            D_fun = casadi.Function('D_fun',{var_fun},{D});
            D = D_fun([x;u]);
        end

        function G = G(obj,x,u) % derivative wrt w
            %TODO, it would need to be PD!
        end
        % add constraints

    end

end

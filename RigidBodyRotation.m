% File: RigidBodyRotation.m
% Author: Antoine Leeman
% Date: 02 Feb 2023
% License: MIT License
% Link: https://arxiv.org/abs/2301.04943
% Citation:
% @article{leeman2023robust,
%   title={Robust Nonlinear Optimal Control via System Level Synthesis},
%   author={Leeman, Antoine P and K{\"o}hler, Johannes and Zanelli, Andrea and Bennani, Samir and Zeilinger, Melanie N},
%   journal={arXiv preprint arXiv:2301.04943},
%   year={2023}
% }
% Description: Class of the dynamical system considered.

classdef RigidBodyRotation < DynamicalSystem
    
    properties

    end
    
    methods
            % constructor function
        function obj = RigidBodyRotation()
            obj@DynamicalSystem(7,3,3,0.1)
            obj.parameters.mu = [3.6993    3.7035    3.7170    3.6350    0.6492    4.6086    5.6358];
            obj.parameters.Inertia = diag([5,2,1]);
            obj.linear = false;
        end        
        
        function dt = ode(obj,x,u) % equation of motion of the dynamical system in continuous time (x_dot = ode(x,u) )
            I = obj.parameters.Inertia;
            q = x(1:4);
            w = x(5:7);
            dt = [0.5 * quat_mut([0 w']',q);
                I\(u - cross(w,I*w))];
        end
        
        function [max_mu] = compute_mu(obj,n_points)
            % estimation of the value of mu, via sampling
            rng(0,'twister');
            tic
            M = [ones(1,4), ones(1,3)*obj.w_max, ones(1,3)*obj.T_max];
            max_mu = zeros(1,obj.nx);
            parfor i = 1:n_points %assume symmetrical constraints
                eval = M.*(2*rand(1,10)-1);
                eval(1:4) = eval(1:4)/norm(eval(1:4)); % (!) non-uniform sampling
                max_mu = max(max_mu,obj.eval_mu(eval));
            end
            toc;
            
        end

        function mu = eval_mu(obj,xu)
            import casadi.*
            x_fun = SX.sym('x',obj.nx);
            u_fun = SX.sym('u',obj.nu);
            var_fun = [x_fun;u_fun];
            H = jacobian(jacobian(obj.ddyn(x_fun,u_fun,'rk4'), var_fun), var_fun);
            H_fun = casadi.Function('H_fun',{var_fun},{H});
            
            H = permute(reshape(full(H_fun(xu)),[obj.nx,obj.nx+obj.nu,obj.nx+obj.nu]),[3,2,1]);
            d = size(H);
            for i = 1:d(3)
                mu(i) = 0.5*sum(sum(abs(H(:,:,i))));
            end
        end
        
    end
end

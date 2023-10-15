classdef KKT_SLS < OCP
    %All the KKT related functions for the SLS problem
    
    properties
        current_bo;
        current_adj_corr;
        solver_forward;
    end
    
    methods
        function obj = KKT_SLS(N,Q,R,m,x0,Qf)
            obj@OCP(N,Q,R,m,x0,Qf);
            obj.current_bo = zeros(m.ni,N);
            obj.current_adj_corr = 0;
            obj.initialize_solver_forward('qpoases');
        end
        
        function initialize_solver_forward(obj,solver)
            import casadi.*
            m=obj.m;
            x = casadi.SX.sym('x',m.nx);
            u = casadi.SX.sym('u',m.nu);
            S_cons = [m.A(x,u), m.B(x,u), eye(m.nx) ;...
                m.C(x,u), m.D(x,u), zeros(m.ni,m.nx)];

            % S_cons_ = blkdiag(eye(m.nx),zeros(m.ni,m.nu) );
            % Z_block = kron(diag(ones(1,obj.N),-1), eye(m.nx+m.nu));
            % 
            % A_mat = blkdiag(kron(eye(obj.N), S_cons),zeros(m.nx)) + kron(eye(obj.N), S_cons_)*Z_block';% WRONG!!
            % A_mat = DM(sparse(A_mat));
            % 
            A_mat = zeros(0,m.nx);
            for kk=1:obj.N-1
                A_mat = [[A_mat, zeros((kk-1)*(m.nx+m.ni), m.nu + m.nx)]; zeros(m.nx+m.ni, (kk-1)*(m.nx+m.nu)), S_cons ];
            end
            A_mat = sparsify(DM(A_mat));

            S_cost = blkdiag(obj.Q, obj.R); % todo: add terminal cost

            H_mat = blkdiag(kron(eye(obj.N-1), S_cost),obj.Qf);
            H_mat = DM(sparse(H_mat));

            options =struct;
            obj.solver_forward = conic('solver', solver, struct('a',A_mat.sparsity(), 'h',H_mat.sparsity()),options);            
        end
        
        function [x_bar, u_bar, lambda, mu] = forward_solve(obj)
            
            outputArg = obj.Property1 + inputArg;
        end
    end
end


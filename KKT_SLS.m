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
        
        function initialize_solver_forward(solver)
            m=obj.m;
            S_cons = [m.A_S(), m.B_S(), sparsify(eye(m.nx)),...
                m.C_S(), m.D_S, sparsify(zeros(m.nx))];
            A_mat = kron(m.N, S_cons);

            S_cost = blkdiag(obj.Q, obj.R); % todo: add terminal cost
            H_mat = kron(m.N, S_cost);

            options =struct;
            obj.solver_forward = conic('solver', solver, struct('a',A_mat.sparsity(), 'h', H_mat.sparsity()),options);            
        end
        
        function [x_bar, u_bar, lambda, mu] = forward_solve(obj)
            
            outputArg = obj.Property1 + inputArg;
        end
    end
end


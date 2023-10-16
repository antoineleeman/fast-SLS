classdef KKT_SLS < OCP
    %All the KKT related functions for the SLS problem
    
    properties
        current_bo;
        initial_bo; % without constraint tightening, i.e., nominal
        current_adj_corr;
        solver_forward;
        current_nominal;
        A_mat_fun;
        H_mat; % could be made a function, for instance for non-qaudratic cost
        H_mat_sparsity;
        lbg;
        lbx_fun;
        ubx_fun;
    end
    
    methods
        function obj = KKT_SLS(N,Q,R,m,x0,Qf)
            obj@OCP(N,Q,R,m,x0,Qf);
            obj.current_bo = zeros(m.ni,N);
            obj.current_adj_corr = 0;
            obj = obj.initialize_solver_forward('qpoases');

            import casadi.*
            x = casadi.SX.sym('x',m.nx);
            u = casadi.SX.sym('u',m.nu);

            [~,f] = m.cons(x,u);
            [~,f_f] = m.cons_f(x);
            
            obj.initial_bo = [kron(ones(obj.N,1), [zeros(m.nx,1); f])]; %assume time invariant constraints + move to initialization
            obj.current_bo = obj.initial_bo;

        end
        
        function obj = initialize_solver_forward(obj,solver)
            import casadi.*
            m=obj.m;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            ni_x = m.ni_x;
            N = obj.N;
            y = casadi.SX.sym('y',(N+1)*nx+(N)*nu);
            jj=1;

            A_mat = zeros(0,m.nx);
            for kk=1:obj.N % no constraint on last time step x_N
                var = y(jj:jj+nx+nu-1 );
                jj = jj+nx+nu;
                x = var(1:nx);
                u = var(nx+1:end);
                S = obj.blockConstraint(x,u);
                S_fun = casadi.Function('S',{var},{S});
                columnPadding = casadi.SX.zeros((kk-1)*(m.nx+m.ni), m.nu + m.nx);
                rowPadding = casadi.SX.zeros(m.nx+m.ni, (kk-1)*(m.nx+m.nu));
                A_mat = sparsify(vertcat(horzcat(A_mat, columnPadding), ...
                    horzcat(rowPadding, S_fun(var))));
            end
            obj.A_mat_fun = casadi.Function('A_mat_fun',{y},{A_mat});

            S_cost = blkdiag(obj.Q, obj.R);
            obj.H_mat = blkdiag(kron(eye(obj.N), S_cost),obj.Qf);
            obj.H_mat = DM(sparse(obj.H_mat));

            obj.H_mat_sparsity = obj.H_mat.sparsity();
            options = struct;
            obj.solver_forward = conic('solver', solver, struct('a',A_mat.sparsity(), 'h',obj.H_mat_sparsity),options);
            %obj.lbg = -casadi.SX.inf((N)*(ni+nx),1);
            obj.lbg = [kron(ones(obj.N,1), [zeros(m.nx,1);  -casadi.DM.inf(ni,1)])];
            obj.current_nominal = zeros(N*(nx+nu)+nx,1);

            x0 = casadi.SX.sym('x0',nx);
            ubx = [x0; casadi.DM.inf(N*(nx+nu),1)];
            lbx = [x0; -casadi.DM.inf(N*(nx+nu),1)];

            obj.ubx_fun= casadi.Function('ubx_fun',{x0},{ubx});
            obj.lbx_fun= casadi.Function('lbx_fun',{x0},{lbx});
        end

        function [x_bar, u_bar, lambda_bar, mu_bar] = forward_solve(obj,x0)
            import casadi.*
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;

            A_current = obj.A_mat_fun(obj.current_nominal);
            sol = obj.solver_forward('a',A_current,'h',obj.H_mat,'lba',obj.lbg,'uba',obj.current_bo,'g',obj.current_adj_corr ,'lbx',obj.lbx_fun(x0),'ubx',obj.ubx_fun(x0));
            
            y_sol = reshape([sol.x;zeros(nu,1)], [nx+nu, N+1]);
            x_bar = y_sol(1:nx,:);
            u_bar = ysol(nx+1:end, 1:N);
            %does it speed up to add constraints on x and u?
            % need to add the initial conditions!
        end

        function S_cons = blockConstraint(obj,x,u)
            import casadi.*
            m=obj.m;
            S_cons = [m.A(x,u), m.B(x,u), -eye(m.nx) ;...
                m.C(x,u), m.D(x,u), zeros(m.ni,m.nx)];
        end
    end
end


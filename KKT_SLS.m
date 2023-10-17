classdef KKT_SLS < OCP
    %All the KKT related functions for the SLS problem
    
    properties
        current_bo;
        current_adj_corr;
        solver_forward;
        current_nominal;

        dh_dyk;
        dh_dxN;
        H_mat; % could be made a function, for instance for non-qaudratic cost
        H_mat_sparsity;
        lbg;
        ubg_current;
        nominal_ubg;
        lbx_fun;
        ubx_fun;

        Q_reg;
        R_reg;

        beta_kj;
        epsilon;
        eta_kj;
        mu_current;
        K_current;

        A_dyn_current;
        B_dyn_current;
        C_cons_current;
        D_cons_current;
    end
    
    methods
        function obj = KKT_SLS(N,Q,R,m,x0,Qf)
            obj@OCP(N,Q,R,m,x0,Qf);
            obj.current_bo = zeros(m.ni,N+1); % there should not be a bo for the last input!
            obj.current_adj_corr = 0;
            obj.epsilon = 1e-3;
            obj = obj.initialize_solver_forward('qpoases');% try osqp!
            obj.Q_reg = 1e-3*eye(m.nx);
            obj.R_reg = 1e-3*eye(m.nu);
            
            import casadi.*
            x = casadi.SX.sym('x',m.nx);
            u = casadi.SX.sym('u',m.nu);

            [~,f] = m.cons(x,u);
            [~,f_f] = m.cons_f(x);
            
            obj.nominal_ubg = [kron(ones(obj.N,1), [zeros(m.nx,1); f])]; % assume time invariant constraints + move to initialization
            obj.ubg_current = obj.nominal_ubg;

        end
        
        function obj = initialize_solver_forward(obj,solver)
            import casadi.*
            m=obj.m;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            %ni_x = m.ni_x;
            N = obj.N;
            y = casadi.SX.sym('y',(N+1)*nx+(N)*nu);
            jj=1;

            A_mat = zeros(0,nx);

            obj.A_dyn_current = cell(1,N);
            obj.B_dyn_current = cell(1,N);
            obj.C_cons_current = cell(1,N);
            obj.D_cons_current = cell(1,N);

            obj.current_nominal = zeros(N*(nx+nu)+nx,1);
            current_nominal_mat = reshape([sol.x;zeros(nu,1)], [nx+nu, N+1]);

            for kk=1:obj.N % no constraint on last time step x_N: could add a equality or box constraint
                var = y(jj:jj+nx+nu-1 );
                jj = jj+nx+nu;
                xkk = current_nominal_mat(1:nx,kk);
                ukk = current_nominal_mat(nx+1:end,kk);

                obj.A_dyn_current{kk} = m.A(xkk,ukk);
                obj.B_dyn_current{kk} = m.B(xkk,ukk);
                obj.C_cons_current{kk} = m.C(xkk,ukk);
                obj.D_cons_current{kk} = m.D(xkk,ukk);
                
                S_fun =[obj.A_dyn_current{kk}, obj.B_dyn_current{kk}, -eye(m.nx) ;...
                obj.C_cons_current{kk}, obj.D_cons_current{kk}, zeros(m.ni,m.nx)];
                columnPadding = casadi.SX.zeros((kk-1)*(m.nx+m.ni), m.nu + m.nx);
                rowPadding = casadi.SX.zeros(m.nx+m.ni, (kk-1)*(m.nx+m.nu));
                
                A_mat = sparsify(vertcat(horzcat(A_mat, columnPadding), ...
                    horzcat(rowPadding, S_fun)));
                % affine part in the dynamics?
            end

            S_cost = blkdiag(obj.Q, obj.R);
            obj.H_mat = blkdiag(kron(eye(obj.N), S_cost),obj.Qf);
            obj.H_mat = DM(sparse(obj.H_mat));

            obj.H_mat_sparsity = obj.H_mat.sparsity();
            options = struct;
            obj.solver_forward = conic('solver', solver, struct('a',A_mat.sparsity(), 'h',obj.H_mat_sparsity),options);
            obj.lbg = [kron(ones(obj.N,1), [zeros(m.nx,1);  -casadi.DM.inf(ni,1)])];

            x0 = casadi.SX.sym('x0',nx);
            ubx = [x0; casadi.DM.inf(N*(nx+nu),1)];
            lbx = [x0; -casadi.DM.inf(N*(nx+nu),1)];

            obj.ubx_fun= casadi.Function('ubx_fun',{x0},{ubx});
            obj.lbx_fun= casadi.Function('lbx_fun',{x0},{lbx});

            % initialization of beta
            obj.beta_kj = obj.epsilon.*ones(ni,N);
            obj.eta_kj = obj.epsilon.*casadi.DM.zeros(ni,N);

        end

        function [obj, x_bar, u_bar, lambda_bar, mu_bar] = forward_solve(obj,x0)
            import casadi.*
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;

            sol = obj.solver_forward('a',obj.A_current,'h',obj.H_mat,'lba',obj.lbg,'uba',obj.ubg_current,'g',obj.current_adj_corr ,'lbx',obj.lbx_fun(x0),'ubx',obj.ubx_fun(x0));
            obj.A_current = obj.A_mat_fun(sol.x);
            obj.current_nominal = sol.x;

            y_sol = reshape([sol.x;zeros(nu,1)], [nx+nu, N+1]);
            x_bar = y_sol(1:nx,:);
            u_bar = y_sol(nx+1:end, 1:N);

            dual = reshape(sol.lam_a, [nx+ni,N]);            
            mu_bar = dual(nx+1:end,:);
            lambda_bar = [sol.lam_x(1:nx), dual(1:nx,:)];

            obj.mu_current = mu_bar;
        end

        function [obj, K] = backward_solve(obj) % could be computed using HPIPM! + we are evaluating twice the dynamics
            import casadi.*
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            
            S = cell(N+1,1);
            K = cell(N,1);

            %make is a generic function in the initialization, so that
            %there is no need to redefine the functions each time.
            import casadi.*
            x = casadi.SX.sym('x',m.nx);
            u = casadi.SX.sym('u',m.nu);
            var = [x;u];
            C_f = m.Cf(x);
            C_f_fun = casadi.Function('C_f_fun',{x},{C_f});
            

            C = m.C(x,u);
            C_fun = casadi.Function('C_fun',{var},{C});
            
            D = m.D(x,u);
            D_fun = casadi.Function('D_fun',{var},{D});

            y_current = reshape([obj.current_nominal;zeros(nu,1)], [nx+nu, N+1]);
            xf = y_current(1:nx,end);

            S{N+1} = C_f_fun(xf)' * diag(obj.eta_kj(1:m.ni_x,N)) *C_f_fun(xf)+ obj.Q_reg;

            for kk=N:-1:1
                ykk = y_current(:,kk);
                xkk = ykk(1:nx);
                ukk = ykk(1:nu);

                C = C_fun(ykk);
                D = D_fun(ykk);

                Ck = [C'; D']*diag(obj.eta_kj(:,kk))*[C, D];

                Cxk = Ck(1:nx, 1:nx) + obj.Q_reg;
                Cuk = Ck(nx+1:end, nx+1:end) + obj.R_reg;
                % Cxuk missing!
                A = m.A(xkk,ukk);
                B = m.B(xkk,ukk);
                K{kk} = -(Cuk+B'*S{kk+1}*B)\(B'*S{kk+1}*A);
                S{kk} = Cxk + A'*S{kk+1}*A + A'*S{kk+1}*B*K{kk};
            end
            obj.current_K = K;

            
        end
        % 
        % function [obj, beta] = compute_backoff(obj)
        %     import casadi.*
        %     m=obj.m;
        %     N = obj.N;
        %     nx = m.nx;
        %     nu = m.nu;
        %     ni = m.ni;
        %     Phi_kj = cell(N+1,1);
        %     Phi_kj{1} = eye(nx);
        %     for kk=1:N
        %         Phi_kj{kk+1} = m.A
        % 
        %     end
        % end

        function S_cons = blockConstraint(obj,x,u)
            import casadi.*
            m=obj.m;
            S_cons = [m.A(x,u), m.B(x,u), -eye(m.nx) ;...
                m.C(x,u), m.D(x,u), zeros(m.ni,m.nx)];
        end

        function obj = solve_beta_stationnarity(obj)
            for kk=1:obj.N
                obj.eta_kj(:,kk) = obj.mu_current(kk)./sqrt(obj.beta_kj(:,kk));
            end
        end


    end
end


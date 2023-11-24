%udpate !!


classdef inexact_KKT_SLS < KKT_SLS
    %All the KKT related functions for the SLS problem
    % re-implementation so that the Riccati recursions are solved inexactly
    
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
        A_current;

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
    end
    
    methods
        function obj = inexact_KKT_SLS(N,Q,R,m,Qf)
            obj@OCP(N,Q,R,m,Qf);
            obj.current_bo = zeros(m.ni,N+1); % there should not be a bo for the last input!
            obj.current_adj_corr = 0;
            obj.epsilon = 1e-3;
           %obj = obj.initialize_solver_forward('qpoases');% try osqp! + HPIPM
            obj = obj.initialize_solver_forward('osqp');
            %obj = obj.initialize_solver_forward('ipopt');

            obj.Q_reg = 1e-3*eye(m.nx);
            obj.R_reg = 1e-3*eye(m.nu);
            
            % import casadi.*
            % x = casadi.SX.sym('x',m.nx);
            % u = casadi.SX.sym('u',m.nu);
            
            obj.nominal_ubg = [kron(ones(obj.N,1), [zeros(m.nx,1); m.d])]; % assume time invariant constraints + move to initialization
            obj.ubg_current = obj.nominal_ubg;

        end
        
        function [feasible,ii, x_bar, u_bar, K] = solve(obj,x0) %% initial conditions should be given here: after the initialization
            %parrallelized, or not?
            %inexact or not?
            % catch HPIPM errors
            % check size of x0
            
            MAX_ITER = 10;
            CONV_EPS = 1e-6;
            m = obj.m;
            N = obj.N;
            current_x = zeros(m.nx,N+1);
            current_u = zeros(m.nu,N);


            % figure(1);
            % clf;
            % hold on;
            % rectangle('Position',[-5 -5 10 10]); %constraints
            % axis equal;
            % 
            % figure(2);
            % clf;
            it_x = cell(MAX_ITER,1);
            it_u = cell(MAX_ITER,1);

            feasible = false;
            for ii=1:MAX_ITER

                try 
                    [obj, x_bar, u_bar, lambda_bar, mu_bar] = obj.forward_solve(x0);
                catch e
                    if contains(e.message, 'error_on_fail')
                        disp('infeasible forward solve -- use soft constraints');
                        return;
                    else
                        rethrow(e);
                    end
                end
                
                it_x{ii} = x_bar;
                it_u{ii} = u_bar;
                if full(max(max(max(current_x-x_bar)),max(max(current_u-u_bar)))) <= CONV_EPS
                    disp('converged!');
                    feasible =true;
                    % if ii>4
                    %     disp('long convergence');
                    %     for ii = 1:9
                    %         tab(ii) = norm(full(it_x{ii} -it_x{ii+1}) ,'inf')
                    %     end
                    % end
                    return;
                else
                    current_x = x_bar;
                    current_u = u_bar;
                end
                obj = obj.update_cost_tube();
                [obj, K] = obj.backward_solve();


                [obj,beta] = obj.update_backoff();
                %[obj,~] = obj.update_backoff_parr();
            end
            disp('no feasible solution found');

        end


        function obj = initialize_solver_forward(obj,solver)
            import casadi.*
            m=obj.m;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            %ni_x = m.ni_x;
            N = obj.N;
            obj.A_dyn_current = cell(1,N);
            obj.B_dyn_current = cell(1,N);
            obj.C_cons_current = cell(1,N);

            obj.current_nominal = zeros(N*(nx+nu)+nx,1);
            current_nominal_mat = reshape([obj.current_nominal;zeros(nu,1)], [nx+nu, N+1]);
            A_mat = zeros(0,nx);

            for kk=1:obj.N % no constraint on last time step x_N: could add a equality or box constraint
                xkk = current_nominal_mat(1:nx,kk);
                ukk = current_nominal_mat(nx+1:end,kk);

                % we should not call this function multiple times, but
                % rather save it for later use.
                obj.A_dyn_current{kk} = m.A;
                obj.B_dyn_current{kk} = m.B;
                obj.C_cons_current{kk} = m.C;
                
                S_fun =[obj.A_dyn_current{kk}, obj.B_dyn_current{kk}, -eye(m.nx) ;...
                obj.C_cons_current{kk}, zeros(m.ni,m.nx)];
                columnPadding = casadi.DM.zeros((kk-1)*(m.nx+m.ni), m.nu + m.nx);
                rowPadding = casadi.DM.zeros(m.nx+m.ni, (kk-1)*(m.nx+m.nu));
                
                A_mat = sparsify(vertcat(horzcat(A_mat, columnPadding), ...
                    horzcat(rowPadding, S_fun)));
                % affine part in the dynamics?
            end
            obj.A_current = A_mat;

            S_cost = blkdiag(obj.Q, obj.R);
            obj.H_mat = blkdiag(kron(eye(obj.N), S_cost),obj.Qf);
            obj.H_mat = DM(sparse(obj.H_mat));

            obj.H_mat_sparsity = obj.H_mat.sparsity();
            options = struct;
            obj.solver_forward = conic('solver', solver, struct('a',A_mat.sparsity(), 'h',obj.H_mat_sparsity),options);
            obj.lbg = [kron(ones(obj.N,1), [zeros(m.nx,1);  -casadi.DM.inf(ni,1)])];

            x0 = casadi.SX.sym('x0',nx);
            ubx = [x0+obj.epsilon; casadi.DM.inf(N*(nx+nu),1)];
            lbx = [x0-obj.epsilon; -casadi.DM.inf(N*(nx+nu),1)];

            obj.ubx_fun= casadi.Function('ubx_fun',{x0},{ubx});
            obj.lbx_fun= casadi.Function('lbx_fun',{x0},{lbx});

            % initialization of beta
            obj.beta_kj = obj.epsilon.*ones(ni,N); % need all the other columns ..
            obj.eta_kj = obj.epsilon.*casadi.DM.zeros(ni,N);

        end

        function [obj, x_bar, u_bar, lambda_bar, mu_bar] = forward_solve(obj,x0)
            import casadi.*
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;

            c = obj.ubg_current;

            sol = obj.solver_forward('a',obj.A_current,'h',obj.H_mat,'lba',obj.lbg,'uba',obj.ubg_current,'g',obj.current_adj_corr ,'lbx',obj.lbx_fun(x0),'ubx',obj.ubx_fun(x0));


            %obj.lbx_fun(x0) and should also include the terminal condition
            
            
            %obj.A_current = obj.A_mat_fun(sol.x);
            
            % add update A_current;
            obj.current_nominal = sol.x;

            y_sol = reshape([sol.x;zeros(nu,1)], [nx+nu, N+1]);%make this as a function
            x_bar = y_sol(1:nx,:);
            u_bar = y_sol(nx+1:end, 1:N);

            dual = reshape(sol.lam_a, [nx+ni,N]);            
            mu_bar = dual(nx+1:end,:);
            lambda_bar = [sol.lam_x(1:nx), dual(1:nx,:)];

            obj.mu_current = full(mu_bar);
        end

        function [obj, K] = backward_solve(obj) % could be computed using HPIPM! + we are evaluating twice the dynamics

            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            
            S = cell(N+1,1);
            K = cell(N,1);

            C_f = m.Cf;
            C = m.C;

            y_current = reshape([obj.current_nominal;zeros(nu,1)], [nx+nu, N+1]);
            [x_bar, u_bar] = obj.convert_y_to_xu(m,obj.current_nominal);
            
            xf = x_bar(1:nx,end);

            S{N+1} = C_f' * diag(obj.eta_kj(1:m.ni_x,N)) *C_f+ obj.Q_reg;

            for kk=N:-1:1
                xkk = x_bar(kk);
                ukk = u_bar(kk);

                %C = obj.C_cons_current{kk};
                Ck = C'*diag(obj.eta_kj(:,kk))*C;

                Cxk = Ck(1:nx, 1:nx) + obj.Q_reg;
                Cuk = Ck(nx+1:end, nx+1:end) + obj.R_reg;
                
                % Cxuk missing!
                A = obj.A_dyn_current{kk};
                B = obj.B_dyn_current{kk};
                K{kk} = -(Cuk+B'*S{kk+1}*B)\(B'*S{kk+1}*A);
                S{kk} = Cxk + A'*S{kk+1}*A + A'*S{kk+1}*B*K{kk};
            end
            obj.K_current = K;            
        end
        % 
        function [obj, beta_kj] = udpate_backoff_inexact(obj)
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            Phi_x_kj = cell(N,1);
            Phi_u_kj = cell(N,1); % !! should not be the same size!

            Phi_x_kj{1} = m.E;

            beta_kj = zeros(m.ni,N-1); % no tightening for the first time step
            C = m.C;
            for kk=1:N-1
                Phi_u_kj{kk} = obj.K_current{kk}*Phi_x_kj{kk};

                beta_kj(:,kk) = vecnorm(C*full([Phi_x_kj{kk};Phi_u_kj{kk}]),2,2);
                bo_k(:,kk) = beta_kj(:,kk); % needs an extra loop over j
                Phi_x_kj{kk+1} = (obj.A_dyn_current{kk} + obj.B_dyn_current{kk}*obj.K_current{kk})*Phi_x_kj{kk};
                % calculate beta here + bo here
            end
            Cf = m.Cf;

            % we don't do the tightening on the last stage
            %beta_kj(1:4,N+1) = [vecnorm(Cf*full(Phi_x_kj{N+1}),2,2)]; % HARDCODED VALUE!!! We should have a separate variable beta for the last state. It could have more constraints.
            obj.beta_kj = beta_kj;

            beta_0j = zeros(ni,1);% no tughtening for the first time step
            %update_ubg = reshape([beta_0j,[zeros(nx,N-1); m.d - beta_kj]],[(N)*(ni+nx),1]);
            update_ubg = reshape([zeros(nx,N); m.d - [beta_0j, beta_kj]],[(N)*(ni+nx),1]);

%            new_ubg = obj.nominal_ubg - update_ubg;
            obj.ubg_current = update_ubg;
        end

        function [obj, beta_kj] = update_backoff(obj) %% use the same controller
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;

            C = m.C;
            Phi_x_kj = cell(N,N);
            Phi_u_kj = cell(N,N); % !! should not be the same size!

            beta_kj = zeros(m.ni,N-1); % no tightening for the first time step
            for jj=1:N-1 % iteration on the disturbance number
                Phi_x_kj{jj,jj} = m.E; % could be time-varying (nonlinear case)
                for kk=jj:N-1
                    Phi_u_kj{kk,jj} = obj.K_current{kk}*Phi_x_kj{kk,jj}; %% kk or jj?

                    beta_kj(:,kk) =beta_kj(:,kk) + vecnorm(C*full([Phi_x_kj{kk,jj};Phi_u_kj{kk,jj}]),2,2);
                    A_cl = obj.A_dyn_current{kk} + obj.B_dyn_current{kk}*obj.K_current{kk};
                    Phi_x_kj{kk+1,jj} = A_cl*Phi_x_kj{kk,jj};
                end
                bo_k(:,jj) = beta_kj(:,jj);
            end
            obj.beta_kj = beta_kj;

            beta_0j = zeros(ni,1);% no tughtening for the first time step
            update_ubg = reshape([zeros(nx,N); m.d - [beta_0j, beta_kj]],[(N)*(ni+nx),1]);
            obj.ubg_current = update_ubg;

        end
        % 
        % 
        % function [obj, beta_kj] = update_backoff_parr(obj) %% use the same controller
        %     m=obj.m;
        %     N = obj.N;
        %     nx = m.nx;
        %     nu = m.nu;
        %     ni = m.ni;
        % 
        %     C = m.C;
        % 
        %     % Preallocate memory for Phi_x_kj and Phi_u_kj
        %     Phi_x_kj = cell(N-1, N-1);
        %     Phi_u_kj = cell(N-1, N-1);
        % 
        %     % Temporary variables for parallel computation, each 'jj' will write to a unique third dimension
        %     temp_beta_kj = zeros(m.ni, N-1, N-1);
        % 
        %     parfor jj = 1:N-1
        %         local_beta_kj = zeros(m.ni, N-1); % Local variable for beta_kj computation
        %         temp_Phi_x_kj = m.E; % Local variable for Phi_x_kj{jj,jj}
        % 
        %         for kk = jj:N-1
        %         %for kk = 1:5
        %             temp_Phi_u_kj = obj.K_current{kk} * temp_Phi_x_kj; % Local computation
        % 
        %             local_beta_kj(:, kk) = vecnorm(C * full([temp_Phi_x_kj; temp_Phi_u_kj]), 2, 2);
        % 
        %             A_cl = obj.A_dyn_current{kk} + obj.B_dyn_current{kk} * obj.K_current{kk};
        %             temp_Phi_x_kj = A_cl * temp_Phi_x_kj; % Update for next iteration
        % 
        %             % Store results in cell arrays
        %             Phi_x_kj{kk, jj} = temp_Phi_x_kj;
        %             Phi_u_kj{kk, jj} = temp_Phi_u_kj;
        %         end
        % 
        %         temp_beta_kj(:, :, jj) = local_beta_kj; % Writing to a unique slice
        %     end
        % 
        %     % Aggregate results from temp_beta_kj
        %     beta_kj = sum(temp_beta_kj, 3);
        % 
        %     % Compute bo_k
        %     bo_k = beta_kj(:, 1:N-1);
        % 
        %     obj.beta_kj = beta_kj;
        % 
        %     beta_0j = zeros(ni,1);% no tughtening for the first time step
        %     update_ubg = reshape([zeros(nx,N); m.d - [beta_0j, beta_kj]],[(N)*(ni+nx),1]);
        %     obj.ubg_current = update_ubg;
        % 
        % end

        % function S_cons = blockConstraint(obj,x,u)
        %     import casadi.*
        %     m=obj.m;
        %     S_cons = [m.A(x,u), m.B(x,u), -eye(m.nx) ;...
        %         m.C(x,u), m.D(x,u), zeros(m.ni,m.nx)];
        % end

        function obj = update_cost_tube(obj)
            for kk=1:obj.N-1
                obj.eta_kj(:,kk) = obj.mu_current(:,kk)./sqrt(obj.beta_kj(:,kk));
            end
        end

        function obj = check_stationnarity(obj)


        end

        function [x,u] = convert_y_to_xu(obj,m,y)
            y_mat = reshape([y;zeros(m.nu,1)], [m.nx+m.nu, obj.N+1]);
            x = y_mat(1:m.nx,:);
            u = y_mat(m.nx+1:end, 1:obj.N);
        end

        % function y = convert_xu_to_y(x,u)
        % end


    end
end

classdef inexact_KKT_SLS < KKT_SLS
    properties
        %no new properties
    end
    
    methods

        function [obj, K] = backward_solve(obj) % only solve the first Riccati

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
                [K{kk},S{kk}] = obj.riccati_step(A,B,Cxk,Cuk,S{kk+1,jj});
            end
            obj.K_current = K;            
        end

        function [obj, beta_kj] = update_backoff(obj) %% using the same controller from the first Riccati
            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;

            C = m.C;
            Phi_x_kj = cell(N,N);
            Phi_u_kj = cell(N,N);

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

        function obj = update_cost_tube(obj) %% we only need to calculate the eta of the first column
            for kk=1:obj.N-1
                obj.eta_kj(:,kk) = obj.mu_current(:,kk)./sqrt(obj.beta_kj(:,kk));
            end
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



    end
end


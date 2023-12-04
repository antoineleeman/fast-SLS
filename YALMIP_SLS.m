classdef YALMIP_SLS < OCP
    %All the KKT related functions for the SLS problem
    
    properties
        Q_reg;
        R_reg;
        yalmip_solve;
        solver;
    end
    
    methods
        function obj = YALMIP_SLS(N,Q,R,m,Qf,solver)
            obj@OCP(N,Q,R,m,Qf);

            obj.Q_reg = 1e-3*eye(m.nx);
            obj.R_reg = 1e-3*eye(m.nu);
            obj = obj.initialize_solve(solver);

        end
        
        function [feasible] = solve(obj,x0) %% initial conditions should be given here: after the initialization
            [ V0, errorcode] = obj.yalmip_solve(x0);
            if errorcode == 0
                feasible = true;
                disp('optimal solution found')
            else
                % Solver reports an infeasible solution or an error
                feasible = false;
                disp('Unfeasible');
            end


        end


        function obj = initialize_solve(obj,solver)
                        
            m = obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            A = m.A;
            B = m.B;

            % Define decision variables
            Z = sdpvar(nx, N + 1, 'full'); % State trajectory variables
            V = sdpvar(nu, N, 'full');     % Input trajectory variables
            X0 = sdpvar(nx, 1, 'full');    % Initial state variable
            
            Phi_x = sdpvar( (N + 1) * nx, (N + 1) * nx, 'full');
            Phi_u = sdpvar( (N + 1) * nu, (N + 1) * nx, 'full');
            
            % Construct the sigma matrix
            sigma_seq = kron(eye(N), m.E);
            Sigma_mat = blkdiag(eye(nx),sigma_seq);
            
            % Define the objective function
            objective = Z(:,N+1)'*obj.Qf*Z(:,N+1);
            for k=1:N
                objective = objective + Z(:,k)'*obj.Q*Z(:,k) + V(:,k)'*obj.R*V(:,k);
            end
            
            %objective = objective + norm([kron(eye(N+1),m.Q_cost)* Phi_x;kron(eye(N+1),m.Q_cost)* Phi_u],'fro')^2;
            objective = objective + norm([Phi_x;Phi_u],'fro')^2; %% change regulizer

            % Initialise the constraints
            constraints = X0 == Z(:,1);
            
            % Add structure constraints for Phi_x and Phi_u
            for k = 1 : N
                constraints = [constraints, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (N + 1 - k)*nx)];
            end
            
            for k = 1 : N
                constraints = [constraints, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (N + 1 - k)*nx)];
            end
            
            % Define block downshift operator
            Z_block = kron(diag(ones(1,N),-1), eye(nx));
            ZA_block = Z_block*blkdiag(kron(eye(N), A), zeros(nx, nx));
            ZB_block = Z_block*blkdiag(kron(eye(N), B), zeros(nx, nu));
            Id = eye((N + 1)*nx);
            
            % Add System Level Parametrisation constraint
            constraints = [constraints, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];
            
            % Add initial state constraint
            constraints = [ constraints, Z(:,1)==X0 ];
            
            % Add state dynamics constraints
            for k=1:N
                constraints = [ constraints, Z(:,k+1)==A*Z(:,k)+B*V(:,k)];
            end
            
            % % state+input constraints
            C = m.C;
            d = m.d;
            nFx = length(d);
            for ii = 1:N
                for jj = 1: nFx
                    f = C(jj,:); b = d(jj);
                    LHS = f*[Z(:,ii);V(:,ii)];
                    for kk = 1:ii-1
                        Phi_ki = [ Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx);
                            Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx)];
                        LHS = LHS + norm(f* Phi_ki, 2);
                    end
                    constraints = [constraints, LHS <= b];
                end
            end
            
            % % terminal constraint 
            % Ft= X_f(:,1:nx);
            % bt = X_f(:,nx+1);
            % nFt = length(bt);
            % for jj = 1:nFt
            %     f = Ft(jj,:); b = bt(jj);
            %     LHS = f*Z(:,N+1);
            %     for kk = 1:N
            %         LHS = LHS + norm(f*Phi_x(N*nx+1:(N+1)*nx,kk*nx+1:(kk+1)*nx),2);
            %     end
            %     constraints = [constraints, LHS <= b];
            % end
            
            % input constraint
            % Fu = m.F_u;
            % bu = m.b_u;
            % nFu = length(bu);
            % for ii = 1:N
            %     for jj = 1: nFu
            %         f = Fu(jj,:); b = bu(jj);
            %         LHS = f*V(:,ii);
            %         for kk = 1:ii-1
            %             LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),2);
            %         end
            %         constraints = [constraints, LHS <= b];
            %     end
            % end
            
            options = sdpsettings('verbose',1,'solver',solver);
            %options = sdpsettings('verbose',1);

            obj.yalmip_solve = optimizer(constraints,objective,options,[X0],V(1));

        end

    end
end
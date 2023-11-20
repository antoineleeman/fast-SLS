classdef YALMIP_SLS < OCP
    %All the KKT related functions for the SLS problem
    
    properties
        Q_reg;
        R_reg;
        yalmip_solver;
    end
    
    methods
        function obj = YALMIP_SLS(N,Q,R,m,x0,Qf)
            obj@OCP(N,Q,R,m,x0,Qf);

            obj.Q_reg = 1e-3*eye(m.nx);
            obj.R_reg = 1e-3*eye(m.nu);
            

        end
        
        function [feasible, x_bar, u_bar, K] = solve(obj,x0) %% initial conditions should be given here: after the initialization
            

        end


        function obj = initialize_solve(obj,solver)
                        
                        m = obj.m;
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;

            % Define decision variables
            Z = sdpvar(nx, N + 1, 'full'); % State trajectory variables
            V = sdpvar(nu, N, 'full');     % Input trajectory variables
            X0 = sdpvar(nx, 1, 'full');    % Initial state variable
            U_L = sdpvar(nu, 1, 'full');   % learned-input bound variable
            
            Phi_x = sdpvar( (N + 1) * nx, (N + 1) * nx, 'full');
            Phi_u = sdpvar( (N + 1) * nu, (N + 1) * nx, 'full');
            
            % Construct the sigma matrix
            sigma_seq = kron(eye(N), m.E);
            Sigma_mat = blkdiag(eye(nx),sigma_seq);
            
            % Define the objective function
            objective = Z(:,N+1)'*m.Q_cost*Z(:,N+1);
            for k=1:N
                objective = objective + Z(:,k)'*m.Q_cost*Z(:,k) + V(:,k)'*m.R_cost*V(:,k);
            end
            
            %objective = objective + norm([kron(eye(N+1),m.Q_cost)* Phi_x;kron(eye(N+1),m.Q_cost)* Phi_u],'fro')^2;
            objective = objective + norm([Phi_x;Phi_u],'fro')^2; %% add regulizer

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
            
            % state constraints
            Fx = m.F_x;
            bx = m.b_x;
            nFx = length(bx);
            for ii = 1:N
                for jj = 1: nFx
                    f = Fx(jj,:); b = bx(jj);
                    LHS = f*Z(:,ii);
                    for kk = 1:ii-1
                        LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx), 2);
                    end
                    constraints = [constraints, LHS <= b];
                end
            end
            
            % terminal constraint 
            Ft= X_f(:,1:nx);
            bt = X_f(:,nx+1);
            nFt = length(bt);
            for jj = 1:nFt
                f = Ft(jj,:); b = bt(jj);
                LHS = f*Z(:,N+1);
                for kk = 1:N
                    LHS = LHS + norm(f*Phi_x(N*nx+1:(N+1)*nx,kk*nx+1:(kk+1)*nx),2);
                end
                constraints = [constraints, LHS <= b];
            end
            
            Fu = m.F_u;
            bu = m.b_u;
            nFu = length(bu);
            for ii = 1:N
                for jj = 1: nFu
                    f = Fu(jj,:); b = bu(jj);
                    LHS = f*V(:,ii);
                    for kk = 1:ii-1
                        LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),2);
                    end
                    constraints = [constraints, LHS <= b];
                end
            end
            
            options = sdpsettings('verbose',0,'solver',solver);
            obj.yalmip_solver = optimizer(constraints,objective,options,X0);

        end

    end
end
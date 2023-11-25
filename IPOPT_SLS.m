classdef IPOPT_SLS < OCP
    properties

        nx;
        nu;
    end
    
    methods
        function obj = IPOPT_SLS(N,Q,R,m,Qf)
            obj@OCP(N,Q,R,m,Qf);
            obj.nx=m.nx;
            obj.nu=m.nu;
            obj.xf = zeros(obj.nx,1);
        end
        
        function R = v3_to_R(obj,v3_Phi_x)
            R = [];
            n=1;
            for k = 1:obj.N
                R_= [];
                for j=1:k
                    R_ = [R_, v3_Phi_x{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    R_ = [R_, sparse(obj.nx,obj.nx)];
                end
                R = [R;R_];
            end
            
        end
        function M = v3_to_M(obj,v3_Phi_u)
            M = [];
            n=1;
            for k = 1:obj.N
                M_ = [];
                for j=1:k
                    M_ = [M_, v3_Phi_u{n}];
                    n = n+1;
                end
                for j=k+1:obj.N
                    M_ = [M_, sparse(obj.nu,obj.nx)];
                end
                M = [M;M_];
            end
        end
        
        function R = v3_to_Sigma(obj,v3_Sigma,d)
            R = [];
            n=1;
            nx = obj.nx;
            for k = 1:obj.N
                R_= [];
                for j=1:k-1
                    R_ = [R_, v3_Sigma{n}]; 
                    n = n+1;
                end
                d_k = d((k-1)*nx+1:k*nx);
                R_ = [R_, diag(d_k)];
                for j=k+1:obj.N
                    R_ = [R_, sparse(nx,nx)];
                end
                R = [R;R_];
            end
            
        end
        
        
        function L = Phi_line_k(obj, k, v3_Phi)
            %k is the number of the line, starting from 1
            L = [];
            last = k*(k+1)/2;
            first = k*(k-1)/2+1;
            for j=first:last
                L = [L, v3_Phi{j}];
            end
        end
        
        function AB = buildNonlinearMap(obj,m,Z,V)
            nx = obj.nx;
            nu = obj.nu;
            N = obj.N;
            
            AB = sparse(nx, N*(nx+nu));
            AB(1:nx,1:nx) = eye(nx);
            
            for i =1:N-1
                A = m.A(Z(:,i+1),V(:,i+1));
                B = m.B(Z(:,i+1),V(:,i+1));
                mat_1 = sparse(nx,(i-1)*nx);
                mat_2 = [-A, eye(nx)];%-A(Z(:,i))
                mat_3 =  sparse(nx,(N-2)*nx + (i-1)*(nu-nx));
                mat_4 = -B;
                mat_5 = sparse(nx, (N-i)*nu );
                AB_ = [mat_1, mat_2, mat_3,mat_4,mat_5];
                AB = [AB;AB_];
            end
        end
        
        function J = getObjectiveNominal(obj,m,Z,V)            
            J =0;
            N = obj.N;
            for i=1:N
                J = J+V(:,i)'*m.R_cost*V(:,i) + (Z(:,i)-m.xf)'*m.Q_cost*(Z(:,i)-m.xf);
            end            
            J = J+ (Z(:,N+1)-m.xf)'*obj.Qf*(Z(:,N+1)-m.xf);
            
        end
        
        function [n,g] = getConstraintsDynamic(obj,m,Z,V)
            
            g = [Z(:,1) - m.x0];
            n = m.nx;
            
            for i=1:obj.N
                g = [g; Z(:,i+1) - m.ddyn(Z(:,i),V(:,i))];
            end
            n = n+ obj.N*m.nx;
        end
        
        function [n,g] = getNominalConstraints(obj,m,Z,V)
            C = m.C;
            c = m.d;
            g = [];
            for i =1:m.ni
                for k=1:obj.N
                    g = [g;C(i,:)*[Z(:,k);V(:,k)] - c(i)];
                end
            end
            n = obj.N * m.ni;
        end
        
        
        function [n,g] = getNonlinearMapConstraints(obj, m,Z,V, Phi_x, Phi_u, Sigma, d)
            N = m.N;
            nx = m.nx;
            nu = m.nu;
            
            R = obj.v3_to_R(Phi_x);
            M = obj.v3_to_M(Phi_u);
            Sigma_with_diag = obj.v3_to_Sigma(Sigma,d);
            
            AB = obj.buildNonlinearMap(m,Z,V);
            
            g_full = [reshape(AB*[R;M] - Sigma_with_diag, [N^2*nx^2,1])];
            [n,g] = obj.reducedMap(m,g_full);
        end
        
        function [n,g] = reducedMap(obj,m,g)
            N = obj.N;
            nx = m.nx;
            v_NZ = find(reshape(kron(tril(ones(N),0),ones(nx)), [N^2*nx^2,1])); % remove all trivial equality constraints
            g = g(v_NZ);
            n = N*(N+1)/2*nx^2;
        end
        
        function [one_norm, n,g_ineq, slack, n_var] = one_norm(obj, v)
            
            import casadi.*
            n_slack = length(v);
            slack = MX.sym('slack',n_slack); % 1-norm always computed on vertical vectors
            one_norm = MX.sym('slack',1);
            
            g_ineq = [ v - slack; -v - slack; sum(slack) - one_norm];
            
            n = 2*n_slack +1;
            n_var = n_slack +1;
        end
        
        function [ one_norm_line, n_ineq,g_ineq, var_slack,n_var_slack ] = vec_one_norm(obj, M)            
            var_slack = [];
            n_var_slack = 0;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            one_norm_line = [];
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = obj.one_norm(M(i,:)');
                one_norm_line = [one_norm_line; one_norm];
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var];
                n_var_slack  = n_var_slack + n_var;
            end
        end
        
        function  [ inf_norm, n_ineq,g_ineq, var_slack,n_var_slack ] = mat_inf_norm(obj, M)
            import casadi.*
            inf_norm = MX.sym('slack_tube',1);
            
            var_slack = [];
            n_var_slack = 1;
            n_ineq = 0;
            g_ineq = [];
            dim = size(M);
            for i = 1 : dim(1)
                [one_norm,n,g, var, n_var] = obj.one_norm(M(i,:)');
                n_ineq = n_ineq+n;
                g_ineq = [g_ineq;g];
                var_slack = [var_slack;var;one_norm];
                n_var_slack  = n_var_slack + n_var;
                
                g_ineq = [g_ineq;one_norm - inf_norm];
                n_ineq = n_ineq + 1;
            end
        end
        
        function [n_ineq,g_ineq, var_cons, n_var_cons] = getLinearConstraints(obj,m,Z,V, v3_Phi_x,v3_Phi_u)
            N = m.N;
            C_x = m.F_x;
            D_x = m.b_x;
         
            C_u = m.F_u;
            D_u = m.b_u;
            
            var_cons = [];
            n_var_cons =0;
            n_ineq = 0;
            g_ineq = [];
            for k=1:N
                Phi_k = Phi_line_k(obj, k, v3_Phi_x);
                CPhi_x = C_x * Phi_k;                                
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, CPhi_x);         
                g_ineq = [g_ineq;g];
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm];
                n_var_cons = n_var_cons + n_var_slack;
                g_ineq = [g_ineq; C_x*Z(:,k+1)+one_norm - D_x];
                
                n_ineq = n_ineq + length(D_x);
            end
            
            for k=1:N
                Phi_k = Phi_line_k(obj, k, v3_Phi_u);
                CPhi_u = C_u * Phi_k;                               
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, CPhi_u);         
                g_ineq = [g_ineq;g];
                n_ineq = n_ineq+n;
                var_cons = [var_cons;var;one_norm];
                n_var_cons = n_var_cons + n_var_slack;
                g_ineq = [g_ineq; C_u*V(:,k+1)+ one_norm- D_u];
                
                n_ineq = n_ineq + length(D_u);
            end
            g_ineq = [g_ineq; C_u*V(:,1) - D_u;C_x*Z(:,1) - D_x; ];
            n_ineq = n_ineq + length(D_x) + length(D_u);
            
        end

        
        function [n_ineq,g_ineq,n_slack, slack] = getConstraintFilter(obj,m,Z,V,v3_Phi_x,v3_Phi_u,d,delta)
            nx = m.nx;
            v_P = vecnorm(m.E,1,2);
            g_ineq = [];
            n_ineq = 0;
            
            slack = [];
            n_slack =0;
               
            for i=1:size(m.theta_v,2)
                k =1 ;
                d_k = d((k-1)*nx+1:k*nx);
                
                theta = m.theta_v(:,i);
                [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, m.ddyn_theta(Z(:,1), V(:,1))*theta);
                g_ineq = [g_ineq; g; one_norm + v_P - d_k];
                n_ineq = n_ineq + n +nx;
                
                slack = [slack;var;one_norm];
                n_slack = n_slack + n_var_slack;
                for k = 2:obj.N
                    d_k = d((k-1)*nx+1:k*nx);               
                    Phix_k = obj.Phi_line_k(k-1, v3_Phi_x);
                    Phiu_k = obj.Phi_line_k(k-1, v3_Phi_u);

                    LHS = [m.A_theta(Z(:,k),V(:,k))*theta*Phix_k + m.B_theta(Z(:,k), V(:,k))*theta*Phiu_k,...
                        m.ddyn_theta(Z(:,k), V(:,k))*theta]; 
                    [one_norm, n,g, var, n_var_slack] = vec_one_norm(obj, LHS);
                    g_ineq = [g_ineq;g];
                    n_ineq = n_ineq+n;

                    slack = [slack;var;one_norm];                
                    n_slack = n_slack + n_var_slack;

                    g_ineq = [g_ineq; one_norm + v_P + delta(k-1)^2*m.mu' - d_k];
                    n_ineq = n_ineq + nx;
                end
            end
        end
        
    
        function [n_ineq,g_ineq, n_var_norm, var_norm] = getConstraintTube(obj,m,delta, v3_Phi_x, v3_Phi_u)
            g_ineq = [];
            n_ineq = 0;
            n_var_norm = 0;
            var_norm = [];
            
            for k = 1:obj.N-1
                Phix_k = obj.Phi_line_k(k, v3_Phi_x);
                Phiu_k = obj.Phi_line_k(k, v3_Phi_u);
                Phi_k = [Phix_k;Phiu_k];
                
                [ inf_norm, n_ineq_norm,g_ineq_norm, var_slack_norm,n_var_slack_norm ] = obj.mat_inf_norm(Phi_k);
                
                n_var_norm = n_var_norm + n_var_slack_norm;
                var_norm = [var_norm;var_slack_norm;inf_norm];
                n_ineq = n_ineq + n_ineq_norm;
                g_ineq = [g_ineq; g_ineq_norm];
                
                g_ineq = [g_ineq; inf_norm - delta(k)];
                n_ineq = n_ineq+1;
            end
        end
       
        
        function [y,n] = getVariablesNominal(obj,Z,V)
            N = obj.N;
            nx = obj.nx;
            nu = obj.nu;
            y = [reshape(Z,[(N+1)*nx,1]); reshape(V,[(N+1)*nu,1])];
            n = (N+1)*(nx+nu);
        end
        
        function [y,n] = getVariablesTube(obj,delta)
            N = obj.N;
            y = reshape(delta,[N-1,1]);
            n = N-1;
        end
        
        function [y,n] = getVariablesResponses(obj,m,Phi_x,Phi_u)
            y = [];
            N = obj.N;
            nu = m.nu;
            nx = m.nx;
            for i = 1:(N+1)*N/2
                y = [y; reshape(Phi_u{i},[nu*nx,1])];
            end
            
            for i = 1:(N+1)*N/2
                y = [y;reshape(Phi_x{i},[nx^2,1])];
            end
            n = (N+1)*N/2 * nu*nx + (N+1)*N/2*nx^2;
        end
        
       
       function [y,n] = getVariablesFilter_onlydiag(obj,m,d)
            y = [];
            N = obj.N;
            nu = m.nu;
            nx = m.nx;
             y = [d];
            n = N*nx;
       end
         
        function Z = vecToNominalState(obj,m,y)
            nx = m.nx;
            N = obj.N;
            Z = reshape(y,[nx,N+1]);
        end
        
        function V = vecToNominalInput(obj,m,y)
            nu = m.nu;
            N = obj.N;
            V = reshape(y,[nu,N+1]);
        end
       
       function [Phi_x,Phi_u] = vecToResponse(obj,m,y)
            nu = m.nu;
            N = obj.N;
            nx = m.nx;
            k = 1;
            for n =1: (N+1)*N/2
                Phi_u{n} = full(reshape(y(k:k+nu*nx-1),[nu,nx]));
                k = k + nu*nx;
            end
            for n =1: (N+1)*N/2
                Phi_x{n} = full(reshape(y(k:k+nx*nx-1),[nx,nx]));
                k = k + nx^2;
            end
       end
        
       function [Sigma, d] = vecToSigma(obj,m,y)
           nu = m.nu;
           N = obj.N;
           nx = m.nx;
           k=1;
           d = y(k : k+nx*N-1);
           Sigma = cellmat(1,N*(N-1)/2,nx,nx);
       end
    end

    
    end
end
classdef par_KKT_SLS < KKT_SLS
        
    methods
    

        function [obj, K] = backward_solve(obj) % could be computed using HPIPM! + we are evaluating twice the dynamics

            m=obj.m;
            N = obj.N;
            nx = m.nx;
            nu = m.nu;
            ni = m.ni;
            
            S = cell(N,N);
            K = cell(N,N);

            C_f = m.Cf;
            C = m.C;

            eta_kj = obj.eta_kj;
            for jj=1:N-1
                eta_Nj = eta_kj{N-1,jj}(1:1:m.ni_x);
                S{N,jj} = C_f' * diag(eta_Nj) *C_f+ obj.Q_reg; %should be Q_f
                A = obj.A_dyn_current{jj};
                B = obj.B_dyn_current{jj};

                for kk=N-1:-1:jj 
                    Ck = C'*diag(obj.eta_kj{kk,jj})*C;
                    Cxk = Ck(1:nx, 1:nx) + obj.Q_reg;
                    Cuk = Ck(nx+1:end, nx+1:end) + obj.R_reg;
                    
                    % Cxuk missing!
                    [K{kk,jj},S{kk,jj}] = obj.riccati_step(A,B,Cxk,Cuk,S{kk+1,jj});
                end
            end
            obj.K_current = K;

        end
    end
end


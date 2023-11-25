classdef RTI_fast_SLS < KKT_SLS

    % no backoff update
    methods
        function [feasible,ii, x_bar, u_bar, K] = solve(obj,x0) %% initial conditions should be given here: after the initialization
            % check size of x0            
            MAX_ITER = 1;

            m = obj.m;
            N = obj.N;
            ni = m.ni;
            current_x = zeros(m.nx,N+1);
            current_u = zeros(m.nu,N);
            it_x = cell(MAX_ITER,1);
            it_u = cell(MAX_ITER,1);

            feasible = false;
            obj.ubg_current = obj.nominal_ubg;
            obj.beta_kj = obj.epsilon.*ones(N-1,N-1,ni);

            for ii=1:MAX_ITER
                try 
                    [obj, ~,x_bar, u_bar, lambda_bar, mu_bar] = obj.forward_solve(x0);
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
                if full(max(max(max(current_x-x_bar)),max(max(current_u-u_bar)))) <= obj.CONV_EPS
                    disp('converged!');
                    feasible =true;
                    return;
                else
                    current_x = x_bar;
                    current_u = u_bar;
                end
                obj = obj.update_cost_tube(); 
                [obj, K] = obj.backward_solve();
%                [obj,beta] = obj.update_backoff();
            end

        end

    end
end


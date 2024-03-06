classdef OCP
    %OCP with quadratic cost; container for all the parameters of the OCP;
    properties
        N;   % Prediction horizon
        Q;   % State cost matrix
        R;   % Input cost matrix
        Qf;  % Final state cost matrix
        m; % model
        xf;

        Q_reg;
        R_reg;
    end
    
    methods
        function obj = OCP(N,Q,R,m,Qf)
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.m = m;
            obj.xf = zeros(m.nx,1);
            obj.Qf = Qf;

            obj.Q_reg = 1e-3*eye(size(Q));
            obj.R_reg = 1e-3*eye(size(R));
        end
        
        function cost = getStageCost(obj,x,u)
            cost=x'*obj.Q*x + u'*obj.R*u;
        end
        function cost = getfinalCost(obj,x)
            cost = (x)'*obj.Qf*(x);
        end
    end
    methods(Static) % todo: implement with square roots!
        function [K,S] = riccati_step(A,B,Cx,Cu,Sk)
            K = -(Cu+B'*Sk*B)\(B'*Sk*A);
            S = Cx + A'*Sk*A + A'*Sk*B*K;
        end
    end
end
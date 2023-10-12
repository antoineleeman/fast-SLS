classdef OCP
    %OCP with quadratic cost; container for all the parameters of the OCP;
    properties
        N;   % Prediction horizon
        Q;   % State cost matrix
        R;   % Input cost matrix
        Qf;  % Final state cost matrix
        m; % model
        x0;
        xf;
    end
    
    methods
        function obj = OCP(N,Q,R,m,x0,Qf)
            obj.N = N;
            obj.Q = Q;
            obj.R = R;
            obj.m = m;
            obj.x0 = x0;
            obj.xf = zeros(m.nx,1);
            obj.Qf = Qf;
            % add assert with size of Q and R
        end
        
        function cost = getStageCost(obj,x,u)
            cost=x'*obj.Q*x + u'*obj.R*u;
        end
        function cost = getfinalCost(obj,x)
            cost = (x)'*obj.Qf*(x);
        end
    end
end


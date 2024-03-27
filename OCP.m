% File: OCP.m
% Author: Antoine Leeman (aleeman(at)ethz(dot)ch)
% Date: 06th March 2024
% License: MIT
% Reference:
%{
@article{leeman2024fast,
  title={Fast System Level Synthesis: Robust Model Predictive Control using Riccati Recursions},
  author={Leeman, Antoine P and K{\"o}hler, Johannes and Messerer, Florian and Lahr, Amon and Diehl, Moritz and Zeilinger, Melanie N},
  journal={arXiv preprint arXiv:2401.13762},
  year={2024}}
%}
% Link: https://arxiv.org/abs/2401.13762
% -----------------------------------------------------------------------------
%%

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

        CONV_EPS;
        obj.CONV_EPS = 1e-8;

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
    methods(Static)
        function [K,S] = riccati_step(A,B,Cx,Cu,Sk) %todo: include Cxu
            K = -(Cu+B'*Sk*B)\(B'*Sk*A);
            S = Cx + A'*Sk*A + A'*Sk*B*K;
        end
        function [K,S] = riccati_step_cholesky(A,B,Cx,Cu,Sk)  % this is slower
            L = chol(Cu + B' * Sk * B, 'lower');
            M = L \ (B' * Sk * A);
            K = -(L' \ M);
            S = Cx + A' * Sk * A + A' * Sk * B * K;
        end

    end

end

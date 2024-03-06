% File: LinearSystem.m
% Author: Antoine Leeman (aleeman@ethz.ch)
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

classdef (Abstract) LinearSystem < DynamicalSystem
    properties
        A; %dyn, ct
        B; %input
        C; % constraint x
        d;
        Cf; % terminal constraint
        df;
        E; % noise matrix
    end

    methods
        % Constructor
        function obj = LinearSystem(A,B,C,E,dt) % initialization of terminal constraint?
            if nargin > 0
                obj.A = A;
                obj.B = B;
                obj.C = C;
                obj.E = E;
                obj.dt =dt;
                [obj.nx, obj.nu] = size(B);
                [~,obj.nw] = size(E);
                [obj.ni,~] = size(C);
                [obj.ni_x,~] = size(Cf);
            end
        end

        function x_p = ddyn(obj,x,u)
            x_p = obj.A*x+obj.B*u;
        end


        %return the same object with soft constraint
        % todo: change the cost!
        function obj = soften(obj)
            C_soft = [obj.C, -eye(obj.ni)];
            B_soft = [obj.B, zeros(obj.nx, obj.ni)];
            
            obj.C = C_soft;
            obj.B = B_soft;
                
            [obj.nx, obj.nu] = size(B_soft);
            [obj.ni,~] = size(C_soft);


        end
    end


end

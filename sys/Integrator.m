% File: Integrator.m
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
classdef Integrator < LinearSystem

    properties
    end
    methods
        function obj = Integrator()
            obj.nx=2;
            obj.nu=1;
            obj.ni = 6;
            obj.ni_x =6;
            obj.nw = 2;
            obj.dt =1;
            obj.E = 0.3*eye(obj.nx);
            obj.A = [1,1;...
                0,1];
            obj.B = [0.5;1];
            obj.C = [[1,0;...
                -1,0;...
                0,1;...
                0,-1;
                zeros(2)],[zeros(4,1);
                1;
                -1]];
            x_max = 5;
            u_max = 3;

            obj.d = [x_max;x_max;x_max;x_max;u_max ;u_max];
            obj.Cf = [1,0;...
                -1,0;...
                0,1;...
                0,-1;...
                0,0;...
                0,0];
            obj.df = [x_max;x_max;x_max;x_max;x_max;x_max];
        end

    end
end


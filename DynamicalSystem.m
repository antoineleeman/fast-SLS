% File: DynamicalSystem.m
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


classdef (Abstract) DynamicalSystem
    properties
        nx; % number of state variables
        nu; % number of input variables
        nw; % size of the disturbance
        ni; % size of the constraints
        ni_x; %number of terminal constraint
        dt; % time step
        parameters;
    end

    methods
        % Constructor
        function obj = DynamicalSystem(nx, nu,nw, dt)
            if nargin > 0
                obj.nx = nx;
                obj.nu = nu;
                obj.nw = nw;
                obj.dt = dt;
                parameters = struct;
            end
        end

    end

end

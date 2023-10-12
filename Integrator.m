% -----------------------------------------------------------------------------
% File: Integrator.m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date: 15th May 2023
% License: MIT
% Reference:
%{
@inproceedings{leeman2023a,
  title = {Predictive safety filter using system level synthesis},
  year = {2023},
  booktitle = {Proceedings of the 5th Annual Learning for Dynamics and Control Conference, PMLR},
  volume = {211},
  author={Leeman, Antoine P. and K{\"o}hler, Johannes and Benanni, Samir and Zeilinger, Melanie N.},
  pages = {1180-1192},
  doi = {10.3929/ethz-b-000615512}
}
%}
% Link: https://arxiv.org/abs/2212.02111
% -----------------------------------------------------------------------------
classdef Integrator < DynamicalSystem

    properties
        nx = 2; % number of state
        nu = 1; % number of input
        ni = 6;
        nw = 2;
        
        
        A = [1,1;0,1];
        B = [0.5;1];
        
        F_x;
        b_x; 
        
        F_u;
        b_u; 
        
        Bw;
        
        N = 10;
        dt;
        
        Q_cost;
        R_cost;
        
        x_max;
        u_max;
        
    end
    
    methods
        function obj = Integrator()
            obj.Bw = 0.3*eye(obj.nx);
            obj.Q_cost = eye(obj.nx);
            obj.R_cost = 100*eye(obj.nu);
            obj.x_max = 5;
            obj.u_max = 3;
            C = [0,0,-1;0,0,1;1,0,0;-1,0,0;0,1,0;0,-1,0];
            c = [obj.u_max ;obj.u_max ;obj.x_max;obj.x_max;obj.x_max;obj.x_max];
            obj.F_x = C(3:end,1:obj.nx);
            obj.b_x = c(3:end);
            
            obj.F_u = C(1:2,obj.nx+1:end);
            obj.b_u = c(1:2);
        end
        
        
    end
end


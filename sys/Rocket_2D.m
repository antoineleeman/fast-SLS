% File: Rocket_2D.m
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
classdef Rocket_2D < LinearSystem

    properties
       L = 1.4
        m = 1.1
        I 
        g = 9.81
        a 
        k_x = 4
        k_th = 5
    end
    methods
        function obj = Rocket_2D()
            obj.I = obj.m * obj.L^2 / 3;  % Correct calculation of I
            obj.a = obj.m * obj.g * obj.L / (2 * obj.I);  % Correct calculation of a


            obj.nx=6;
            obj.nu=2;

            obj.ni = 16;
            obj.ni_x =12;

            obj.nw = 6; % todo: faster code if we implement for nw < nx
            
            dt = 0.075;
            obj.dt =dt;
            sig_w = 0.05;
            k_x = obj.k_x;
            a = obj.a;
            k_th = obj.k_th;
            I = obj.I;
            obj.E = diag([0, sig_w, 0,0,0,sig_w,0]);
            obj.A =[1, dt,0,  0, 0,  0;
                         0,  1, dt*k_x,  0, 0,  0;
                         0,  0,      1, dt, 0,  0;
                         0,  0,   a*dt,  1, 0,  0;
                         0,  0,      0,  0, 1, dt;
                         0,  0,      0,  0, 0,  1];
            A = obj.A;

            obj.B =[        0,      0;
                            0,      0;
                            0,      0;
                     dt*k_th/I,      0;
                            0,      0;
                            0,     dt];
            B = obj.B;

            Hx = kron(eye(size(A, 1)), [1; -1]);
            hx = [15.0; 15.0; 6.0; 6.0; 25.0; 25.0; 8.0; 8.0; 15.0; 0.0; 6.0; 6.0];


            Hu = kron(eye(size(B, 2)), [1; -1]);
            hu = [ones(2, 1) * 15;
                  ones(2, 1) * 8];

            obj.C = blkdiag(Hx, Hu);
            
            obj.d = [hx; hu];
            obj.Cf = Hx;
            obj.df = hx;
        end

    end
end


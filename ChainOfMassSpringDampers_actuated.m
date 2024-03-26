% File: ChainOfMassSpringDampers_actuated.m
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
classdef ChainOfMassSpringDampers_actuated < LinearSystem

    properties
        M; %number of mass
        mass;
        k_constant;
        d_constant;
    end

    methods
        function obj = ChainOfMassSpringDampers_actuated(M)
            obj.mass = 1;
            obj.k_constant=10;
            obj.d_constant=2;
            obj.dt = 0.1;
            
            obj = initialization(obj,M);

            obj.nw = obj.nx;
           
            obj.E = 0.5*eye(obj.nw);

            u_max = 4; %% TYPO FROM PAPER!
            x_max = 4;

            obj.C = [eye(obj.nx),zeros(obj.nx,obj.nu,1);
                -eye(obj.nx),zeros(obj.nx,obj.nu,1);
                zeros(obj.nu,obj.nx) ,eye(obj.nu);
                zeros(obj.nu,obj.nx) ,-eye(obj.nu);
                ];
            obj.d = [x_max*ones(obj.nx,1);x_max*ones(obj.nx,1);u_max*ones(obj.nu,1);u_max*ones(obj.nu,1);];
            obj.Cf = [eye(obj.nx);
                -eye(obj.nx);
                [zeros(2*obj.nu,obj.nx) ]]; % no terminal constraint
            obj.df = [x_max*ones(2*obj.nx,1);u_max*ones(2*obj.nu,1)];
            obj.ni = 2*obj.nx+2*obj.nu;
            obj.ni_x =2*obj.nx+2*obj.nu;
        end
        
        % taken from: https://gitlab.ethz.ch/ics/NMPC-StabilityAnalysis/-/blob/main/LinearExample_ChainOfMassSpringDampers/main.m?ref_type=heads
        function obj = initialization(obj,M)

            nxx=2;%number of local states
            nx=nxx*M;%number of overall states

            m=obj.mass;
            k=obj.k_constant;
            d=obj.d_constant;

            %define cont.-time system matrices \dot{x}=A_c*x+B_c*u
            A_c=zeros(nx);

            if M ==1
                A_c = [0, 1; -k/m, -b/m];
            elseif M==2
                A_c = [ 0,1, 0,0;
                        -2*k/m, -2*d/m, k/m, d/m;
                        0,0,0,1;
                        k/m, d/m, -k/m, -d/m];
            else %M>2
                %i=1:also connected to ground
                A_c(1,:)=[0,1,zeros(1,nx-nxx)];
                A_c(2,:)=1/m*[-k-k,-d-d,k,d,zeros(1,nx-2*nxx)];
                for i=2:M-1
                    A_c((i-1)*nxx+1,:)=[zeros(1,(i-1)*nxx+1),1,zeros(1,nx-i*nxx)];
                    A_c(i*nxx,:)=1/m*[zeros(1,nxx*(i-2)),k,d,-2*k,-2*d,k,d,zeros(1,nx-nxx*(i+1))];
                end
                A_c(M*nxx-1,:)=[zeros(1,M*nxx-1),1];
                A_c(M*nxx,:)=1/m*[zeros(1,nx-2*nxx),k,d,-k,-d]; 
            end


            %actuation on each mass
            B_c=1/m*[kron(eye(M),[0;1])];
            [obj.A,obj.B]=c2d(A_c,B_c,obj.dt);
            obj.nx = nx;
            obj.nu = obj.nx/2;

        end

    end
end


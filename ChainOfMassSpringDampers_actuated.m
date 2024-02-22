% -----------------------------------------------------------------------------
% File: .m
% Author: Antoine Leeman (aleeman@ethz.ch)
% Date:  2023
% License: MIT
% Reference:
%{


}
%}
% Link: 
% -----------------------------------------------------------------------------
classdef ChainOfMassSpringDampers_actuated < LinearSystem

    properties
        M; %number of masses
        mass;
        k_constant;
        d_constant;%mass,spring, damping constant

    end

    methods
        function obj = ChainOfMassSpringDampers_actuated(M)
            obj.mass = 1;
            obj.k_constant=10;
            obj.d_constant=2;
            obj.dt = 0.5; %sampling time
            
            obj = initialization(obj,M);

            obj.nw = obj.nx;
           
            obj.E = 0.5*eye(obj.nw);

            u_max = 0.5;
            x_max = 4;

            %obj.d = u_max;
            obj.C = [eye(obj.nx),zeros(obj.nx,obj.nu,1);
                -eye(obj.nx),zeros(obj.nx,obj.nu,1);
                zeros(obj.nu,obj.nx) ,eye(obj.nu);
                zeros(obj.nu,obj.nx) ,-eye(obj.nu);
                ];
            obj.d = [x_max*ones(obj.nx,1);x_max*ones(obj.nx,1);u_max*ones(obj.nu,1);u_max*ones(obj.nu,1);]; %% wrong??
            obj.Cf = [eye(obj.nx);
                -eye(obj.nx);
                [zeros(2*obj.nu,obj.nx) ]]; % no terminal constraint
            obj.df = [x_max*ones(2*obj.nx,1);u_max*ones(2*obj.nu,1)];
            obj.ni = 2*obj.nx+2*obj.nu;
            obj.ni_x =2*obj.nx+2*obj.nu;


            % %obj.d = u_max;
            % obj.C = [zeros(2,obj.nx) ,ones(2,1)];
            % obj.d = [u_max*ones(2,1)];
            % 
            % obj.Cf = [zeros(2,obj.nx)];
            % obj.df = [u_max*ones(2,1)];
            % obj.ni = 2;
            % obj.ni_x =2;



        end
        
        % taken from: https://gitlab.ethz.ch/ics/NMPC-StabilityAnalysis/-/blob/main/LinearExample_ChainOfMassSpringDampers/main.m?ref_type=heads
        function obj = initialization(obj,M)

            % 1. Define model- Large scale mass-spring damper system
            
            nxx=2;%number of local states
            if M<3
               error('not programmed') 
            end
            nx=nxx*M;%number of overall states

            m=obj.mass;
            k=obj.k_constant;
            d=obj.d_constant;

            %define cont.-time system matrices \dot{x}=A_c*x+B_c*u
            A_c=zeros(nx);

            %i=1:also contected to ground
            A_c(1,:)=[0,1,zeros(1,nx-nxx)];
            A_c(2,:)=1/m*[-k-k,-d-d,k,d,zeros(1,nx-2*nxx)];
            %1<i<M
            for i=2:M-1
                A_c((i-1)*nxx+1,:)=[zeros(1,(i-1)*nxx+1),1,zeros(1,nx-i*nxx)];
                A_c(i*nxx,:)=1/m*[zeros(1,nxx*(i-2)),k,d,-2*k,-2*d,k,d,zeros(1,nx-nxx*(i+1))];
            end
            %i=M, measured output
            A_c(M*nxx-1,:)=[zeros(1,M*nxx-1),1];
            A_c(M*nxx,:)=1/m*[zeros(1,nx-2*nxx),k,d,-k,-d]; 

            %acutation on last mass
            obj.nu=1;
            %B_c=1/m*[zeros(nx-nxx,1);0;1;];
            B_c=1/m*[kron(eye(M),[0;1])];
            [obj.A,obj.B]=c2d(A_c,B_c,obj.dt);
            obj.nx = nx;
            obj.nu = obj.nx/2;

        end

    end
end

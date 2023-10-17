classdef MyCasadiContainer
    properties
        casadi_function  % Property to hold the CasADi function
    end
    
    methods
        function obj = MyCasadiContainer()
            % Constructor
            import casadi.*
            
            % Defining a simple CasADi function
            x = SX.sym('x');
            y = SX.sym('y');
            f = x^2 + y^2;
            
            % Creating a CasADi Function
            obj.casadi_function = Function('f', {x, y}, {f});
        end
        
        function evaluate(obj, x_val, y_val)
            % Method to evaluate the CasAdi function
            res = obj.casadi_function(x_val, y_val);
            disp('Result:');
            disp(full(res));
        end
    end
end

%chatGPT
function ellipse_points = MinkowskiSumEllipsoids(A_matrices)
    % A_matrices: Cell array containing the matrices A_i defining each ellipsoid
    
    % Initialize the combined matrix A
    A = zeros(2, 2);
    
    % % Sum all the matrices A_i
    % for i = 1:length(A_matrices)
    %     if ~isempty(A_matrices{i})
    %         A = A + A_matrices{i}'*A_matrices{i};
    %     else
    %         break
    %     end
    % end
        % Sum all the matrices A_i
    for i = 1:length(A_matrices)
        if ~isempty(A_matrices{i})
            A = A + A_matrices{i};
        else
            break
        end
    end
    A = A*A';
    
    % Eigen decomposition to get the semi-axes of the ellipsoid
    [V, D] = eig(A);
    
    % Define the angle for the ellipse points
    theta = linspace(0, 2*pi, 100);
    
    % Parametric equation of the ellipse
    ellipse_points = (V * sqrt(D)) * [cos(theta); sin(theta)];
    

end

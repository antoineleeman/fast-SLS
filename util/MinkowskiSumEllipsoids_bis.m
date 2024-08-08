%chatGPT
function minkowski_points = MinkowskiSumEllipsoids_bis(A_matrices)
    numPoints = 30;
    % A_matrices: Cell array containing the matrices A_i defining each ellipsoid
    
    % Initialize the combined matrix A
    A = A_matrices{1}'*A_matrices{1}+ 0.0001*eye(2);

    % Eigen decomposition to get the semi-axes of the ellipsoid
    [V, D] = eig(A);

    % Define the angle for the ellipse points
    theta = linspace(0, 2*pi, numPoints); theta(end) = [];

    % Parametric equation of the ellipse
    ellipse_points = (V * sqrt(D)) * [cos(theta); sin(theta)];

    P = polyshape(ellipse_points');

    % Sum all the matrices A_i
    for i = 2:length(A_matrices)
        if ~isempty(A_matrices{i})
            A = A_matrices{i}'*A_matrices{i} + 0.01*eye(2);

            % Eigen decomposition to get the semi-axes of the ellipsoid
            [V, D] = eig(A);

            % Define the angle for the ellipse points
            theta = linspace(0, 2*pi, numPoints); theta(end) = [];

            % Parametric equation of the ellipse
            ellipse_points = (V * sqrt(D)) * [cos(theta); sin(theta)];

            Pb = polyshape(ellipse_points');
            P = minkowskiSum(P,Pb);

            indices = round(linspace(1, length(P.Vertices) ,numPoints));
            P = polyshape(P.Vertices(indices,:));
        else
            break
        end
    end
    minkowski_points = P.Vertices';

    
    

end

function plotZonotope(G, maxGenerators)
    % Function to plot a 2D zonotope with selected generators
    % G: Input matrix where each column is a generator (2 x m)
    % maxGenerators: Maximum number of generators to keep
    
    % Compute the norms of the generators
    norms = vecnorm(G);
    
    % Sort the norms and get the indices of the top generators
    [~, sortedIndices] = sort(norms, 'descend');
    
    % Select the top generators based on maxGenerators
    selectedIndices = sortedIndices(1:min(maxGenerators, size(G, 2)));
    selectedGenerators = G(:, selectedIndices);
    
    % Generate the points of the zonotope
    numGenerators = size(selectedGenerators, 2);
    vertices = zeros(2, 2^numGenerators);
    
    % Create the zonotope vertices by summing combinations of generators
    for i = 1:2^numGenerators
        binaryCombination = dec2bin(i-1, numGenerators) - '0';
        vertices(:, i) = sum(selectedGenerators .* binaryCombination, 2);
    end
    
    % Get the convex hull of the vertices to plot the zonotope
    k = convhull(vertices(1, :), vertices(2, :));
    
    % Plot the zonotope
    figure;
    plot(vertices(1, k), vertices(2, k), 'b-', 'LineWidth', 2);
    hold on;
    plot(vertices(1, :), vertices(2, :), 'ro');
    xlabel('x');
    ylabel('y');
    title('2D Zonotope');
    grid on;
    axis equal;
    hold off;
end

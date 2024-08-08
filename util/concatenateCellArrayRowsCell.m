%ChatGPT generated:
function concatenatedResults = concatenateCellArrayRows(C)
    % Function to concatenate each row of a cell array containing 2x2 doubles
    % C: Input cell array (nRows x nCols)
    % concatenatedResults: Output cell array (nRows x 1) with concatenated rows
    
    % Get the size of the cell array
    [nRows, nCols] = size(C);
    
    % Initialize a new cell array to store the concatenated results
    concatenatedResults = cell(nRows, 1);
    
    % Iterate over each row
    for i = 1:nRows
        concatenatedResults{i} = C{i,:};
    end
end

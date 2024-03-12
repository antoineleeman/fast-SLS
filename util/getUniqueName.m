function [output] = getUniqueName(str)
    % Get the current date and time
    currentDateTime = datetime('now');

    % Convert the date and time to string format
    currentDateTimeString = char(currentDateTime);

    % Replace the ":" and " " characters with "_"
    currentDateTimeString = strrep(currentDateTimeString, ':', '_');
    currentDateTimeString = strrep(currentDateTimeString, ' ', '_');
    
    output = strcat(currentDateTimeString,'__',str,'.mat');
end
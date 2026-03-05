function [convertedCoords] = hw_convertAllCoordinates(VesLocs3DStr)
    convertedCoords = zeros(length(VesLocs3DStr), 2);
    
    for i = 1:length(VesLocs3DStr)
        [convertedX, convertedY] = convertCoordinates(VesLocs3DStr{i});
        convertedCoords(i, :) = [convertedX, convertedY];
    end
end



function [convertedX, convertedY] = convertCoordinates(inputCoords)
    % Remove all spaces from the input string
    inputCoords = strrep(inputCoords, ' ', '');

    % Split the string into X and Y coordinates
    coords = strsplit(inputCoords, ',');
    x = coords{1};
    y = coords{2};
    
    % % Remove the character before the coordinate to keep only the number
    digitX = str2double(x(2:end)); % Take from the second character
    digitY = str2double(y(2:end));
    
    % Determine the sign based on the character
    if x(1) == 'R'
        convertedX = -digitX;
    elseif x(1) == 'L'
        convertedX = digitX;
    end
    
    if y(1) == 'A'
        convertedY = digitY;
    elseif y(1) == 'P'
        convertedY = -digitY;
    end
   
    % fprintf('Converted coordinates: [%d, %d]\n', convertedX, convertedY);
end
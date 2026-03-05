function convertedValue = hw_convertCoordinateValue(inputStr)
    % Remove all spaces from the input string
    inputStr = strrep(inputStr, ' ', '');

    numberStr = inputStr(2:end);
    
    number = str2double(numberStr);
    
    if inputStr(1) == 'A' || inputStr(1) == 'L' || inputStr(1) == 'H'
        convertedValue = number;
    elseif inputStr(1) == 'P' || inputStr(1) == 'R' || inputStr(1) == 'F'
        convertedValue = -number;
    else
        error('Invalid character for coordinate value');
    end
end
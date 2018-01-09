function value = gemLibraryIsInPath()
% This function checks whether the gem and sgem classes are in matlab path.

    try
        if isequal(toStrings(gem('e',32),30), '2.71828182845904523536028747135') && ...
           isequal(toStrings(sgem('e',32),30), '2.71828182845904523536028747135')
            value = true;
        else
            value = false;
        end
    catch me
        value = false;
    end

end


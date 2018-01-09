function value = gemLibraryIsInPath()
% This function checks whether the gem and sgem classes are in matlab path.
 
    persistent gemLiraryPresent;
 
    if isempty(gemLiraryPresent)
        % the first time, we check if the classes exist
        if exist('gem','class') &&  exist('sgem','class')
            gemLiraryPresent = true;
        else
            gemLiraryPresent = false;
        end
    end
    
    value = gemLiraryPresent;
end


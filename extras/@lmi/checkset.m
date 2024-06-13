function [p,d] = checkset(X)
% OBSOLETE USE CHECK

switch nargout
    case 0
        check(X);
    case 1
        p = check(X);
    case 2
        [p,d] = check(X);
end

function [p,d] = checkset(X)
% OBSOLETE USE CHECK

switch nargout
    case 0
        check(lmi(X));
    case 1
        p = check(lmi(X));
    case 2
        [p,d] = check(lmi(X));
end

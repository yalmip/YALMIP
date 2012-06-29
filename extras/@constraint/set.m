function F = set(varargin)

switch nargin
    case 0
        F = lmi;
    case 1
        F = lmi(varargin{1});
    case 2
        F = lmi(varargin{1},varargin{2});
    case 3
        F = lmi(varargin{1},varargin{1},varargin{3});
    case 4
        F = lmi(varargin{1},[],[],1);
        
    otherwise
end

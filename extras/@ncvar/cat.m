function y = cat(varargin)
%CAT (overloaded)

% Author Johan Löfberg 
% $Id: cat.m,v 1.1 2006-08-10 18:00:19 joloef Exp $  

switch varargin{1}
    case 1
        y = vertcat(varargin{2:end});
    case 2
        y = horzcat(varargin{2:end});       
    otherwise
end
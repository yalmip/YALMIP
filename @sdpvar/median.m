function varargout=median(varargin)
%MEDIAN (overloaded)
%
% M = median(x)
%
% MEDIAN is implemented using the overloaded SORT operator.

% Author Johan Löfberg
% $Id: median.m,v 1.3 2007-07-29 17:32:28 joloef Exp $

x = varargin{1};

if nargin > 1 | min(size(x))>1
    error('SDPVAR/MEDIAN only supports simple 1-D median'),
end

switch length(x)
    case 1
        varargout{1} = x;
    case 2
        varargout{1} = sum(x);
    otherwise
        y = sort(x);
        if even(length(x))
            y1 = extsubsref(y,length(x)/2);
            y2 = extsubsref(y,1+length(x)/2);
            varargout{1} = (y1+y2)/2;
        else
            varargout{1} = extsubsref(y,ceil(length(x)/2));
        end
end
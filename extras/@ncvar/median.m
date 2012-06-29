function varargout=median(varargin)
%MEDIAN (overloaded)
%
% M = median(x)
%
% MEDIAN is implemented using the overloaded SORT operator.

% Author Johan Löfberg
% $Id: median.m,v 1.1 2006-08-10 18:00:21 joloef Exp $

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


%         
% 
% switch class(x)
% 
%     case 'double'
%         y = sort(x)
%         
%                 if even(length(x))
%                     y1 = y(length(x)/2);
%                     y2 = y(1+length(x)/2);
%                     varargout{1} = (y1+y2)/2;
%                 else
%                     varargout{1} = y(ceil(length(x)/2));
%                 end
% 
%     case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
% 
%         if nargin > 1 | min(size(varargin{1}))>1
%             error('SDPVAR/MEDIAN only supports simple 1-D median'),
%         end
% 
%         switch length(x)
%             case 1
%                 varargout{1} = x;
%             case 2
%                 varargout{1} = sum(x);
%             otherwise
%                 y = sort(x);
% 
%                 if even(length(x))
%                     y1 = y(length(x)/2);
%                     y2 = y(1+length(x)/2);
%                     varargout{1} = (y1+y2)/2;
%                 else
%                     varargout{1} = y(ceil(length(x)/2));
%                 end
%             otherwise
%                 error('Strange type on first argument in SDPVAR/MEDIAN');
%         end

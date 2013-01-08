function varargout=length(varargin)

% Author Johan Löfberg 
% $Id: length.m,v 1.2 2004-07-19 13:54:36 johanl Exp $   

F = varargin{1};
varargout{1} = length(F.LMIid);

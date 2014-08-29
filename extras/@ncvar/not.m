function varargout = not(varargin)
%NOT (overloaded)
%   
%    y = not(x)       
%    y = ~x    
%
% Short for y = 1-x.
%
% It is assumed that x is a binary variable (either explicitely declared
% using BINVAR, or constrained using BINARY.)
%
%    See also SDPVAR/OR, SDPVAR/AND, BINVAR, BINARY

varargout{1}=1-varargin{1};

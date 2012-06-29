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

% Author Johan Löfberg 
% $Id: not.m,v 1.1 2006-08-10 18:00:21 joloef Exp $   

varargout{1}=1-varargin{1};

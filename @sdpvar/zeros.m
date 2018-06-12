function sys = zeros(varargin)
% ZEROS Creates overloaded SDPVAR zeros
%
 % zeros(N,M,'like',spvar) generates an N-by-M zero matrix
 % which actually is an SDPVAR, which means its elements can be changed to 
 % SDPVARs
 % 
 % Example
 % X = zeros(2,3,'like',sdpvar)
 % y = sdpvar(1);
 % X(1,3) = y;
 % ShouldBeZero = X - [0 0 y;0 0 0]
 
sys = double2sdpvar(zeros(varargin{1:nargin-2}));
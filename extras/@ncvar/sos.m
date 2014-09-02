function X=sos(X,r)
%SOS Declare sum-of-squares structure
%
% F = sos(p,r)
%
% Input
%  p : SDPVAR object
%  r : Desired rank (optional)
% Output
%  F : Constraint object
%
% Example:
%  Typical usage is
%
%   F = sos(p)
%
%  An experimental feature is to search for
%  low rank decompositions. To search for a 
%  decomposition using at most 3 terms, use
%  a second argument
%
%   F = sos(p,3)
%
%  Note that his feature requires the solver LMIRANK.   
%
% See also SOSEX 

if ~issymmetric(X)% X.n>1 | X.m>1
    error('SOS can only be applied to symmetric polynomial matrices');
end
if nargin<2
    r = inf;
end    
X.typeflag = 11;
X.extra.sosid = yalmip('sosid');
X.extra.rank = r;
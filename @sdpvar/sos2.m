function X=sos2(X,weights)
%SOS2 Declare special ordered set of type 2
%
% F = sos(p,w)
%
% Input
%  p : SDPVAR object
%  w : Adjacency weights
% Output
%  F : CONSTRAINT object

% Author Johan Löfberg 
% $Id: sos2.m,v 1.8 2009-10-08 11:11:06 joloef Exp $  

X.typeflag = 50;
if nargin == 1
    X.extra.sosweights = 1:length(X);
else
    X.extra.sosweights = findOutWeights(X,weights)       
end
X = set(X);

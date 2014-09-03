function X=sos1(X,weights)
%SOS1 Declare special ordered set of type 1
%
% F = sos(p,w)
%
% Input
%  p : SDPVAR object
%  w : Priority weights
% Output
%  F : CONSTRAINT object

X.typeflag = 51;
if nargin == 1
    X.extra.sosweights = 1:length(X);
else
    X.extra.sosweights = findOutWeights(X,weights)       
end
X = lmi(X);

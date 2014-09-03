function X=sos2(X,weights,dim)
%SOS2 Declare special ordered set of type 2
%
% F = sos2(p,w,dim)
%
% Input
%  p   : SDPVAR object
%  w   : Adjacency weights (optional)
%  dim : Operatate along dim for matrix argument. Default 1
%
% Output
%  F : CONSTRAINT object
%
% See also sos1

% Normalize arguments
if nargin < 3
    dim = 1;
end
if nargin < 2
    weights = [];
end

% Easy fix for column case
if dim == 2
    X = sos2(X',weights);
    return
end
if dim > 2
    error('DIM should be 1 or 2');    
end

if length(size(X))>2
    error('SOS2 not support for nD arrays');
end

% Matrix case by recursive calls
if min(size(X))>1
    Y = [];
    for i = 1:size(X,1)
        T.type = '()';
        T.subs{1} = i;
        T.subs{2} = ':';       
        Y = [Y;sos2(subsref(X,T),weights)];
    end
    X = Y;
    return
end

if isempty(weights)
    X.typeflag = 50;
    X.extra.sosweights = 1:length(X);
else
    X.typeflag = 50;
    X.extra.sosweights = findOutWeights(X,weights)           
end

X = lmi(X);
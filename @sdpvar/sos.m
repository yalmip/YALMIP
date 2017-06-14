function X=sos(X,r)
%SOS Declare sum-of-squares structure
%
% F = sos(p)
%
% Input
%  p : SDPVAR object
% Output
%  F : Constraint
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

if nargin<2
    r = inf;
end

if any(isinf(getbase(X)))
    error('You have infinite elements in the polynomial');
end
if any(isnan(getbase(X)))
    error('You have NaN elements in the polynomial');
end

if ~is(X,'symmetric')
    % User supplied a vector
    X = reshape(X,prod(size(X)),1);
    Z = [];
    for i = 1:length(X)
        I.type = '()';
        I.subs = {[i]};
        x = subsref(X,I);
        if isnumeric(x)
            if x < 0
                error('You are trying to enforce a negative constant to be SOS!');
            end
        else
            Z = [Z,sos(x)];
        end
    end
    X = Z;
else
    X.typeflag = 11;
    X.extra.sosid = yalmip('sosid');
    X.extra.rank = r;
    X = lmi(X);
end
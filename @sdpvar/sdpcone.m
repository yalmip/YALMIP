function y = sdpcone(varargin)
%SDPCONE Defines a semidefinite constraint X>=0
%
% Input
%    z       : Linear SDPVAR object.
%    h       : Linear scalar SDPVAR object
%
% Example
%
% A standard SDP constraint normally written as [X>=0] is defined with
%    F = sdpcone(X)
%
% The typical use though is when we want to define a large number of
% constraint without overhead. By stacking N matrices of size nxn, sdpcone
% defines N SDP constraints directly, without any analysis or check for
% symmetry etc
%
% X = sdpvar([5 5 5],[5 5 5]);
% F = sdpcone([X{:}]);
%
% Alternatively
% F = sdpcone(X{:});
%
% See also  @SDPVAR/CONE, @SDPVAR/RCONE

try
    y = [varargin{:}];
catch
    error('All matrices must have the same dimension.')
end

[n,m] = size(y);
if n == m
    y = y>=0;
    return
end
if rem(m,n)
    error('This cannot be N stacked nxn matrices.');
end
y.typeflag = 57;
y=lmi(y);
function y = sdpcone(y)
%SDPCONE Low-level operator to define several SDP constraints in vectorized
%form
%
% Input
%    X       : Linear SDPVAR object of size n^2 x N
%
% Example
%
% The typical use is when we want to define a very large number of LMI
% constraints with reduced overhead (no analysis or check for symmetry etc)
%
% This operator is very specialized and low-level, and not normally used...
%
% X = sdpvar(5);Y = sdpvar(5);
% F = sdpcone([X(:) Y(:)]); % Equivalent to [X>=0, Y>=0]
%
% See also  @SDPVAR/CONE, @SDPVAR/RCONE

y.typeflag = 57;
y=lmi(y);
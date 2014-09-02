function  LinearVariables = depends(x)
%DEPENDS Returns indicies to variables used in an SDPVAR object
%
% i = depends(x)
%
% Input
%    x : SDPVAR object
% Output
%    i : DOUBLE

[mt,variabletype] = yalmip('monomtable');

% Simple linear cases
if ~any(variabletype(x.lmi_variables))
   LinearVariables = x.lmi_variables;
else
    LinearVariables = find(any(mt(x.lmi_variables,:),1));    
end
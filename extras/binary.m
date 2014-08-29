function x = binary(x)
%BINARY Overloaded

if isempty(x)
    x = [];
else
    error('BINARY can only be applied to SDPVAR objects or empty doubles');
end
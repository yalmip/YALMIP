function x = integer(x)
%INTEGER Overloaded

if isempty(x)
    x = [];
else
    error('INTEGER can only be applied to SDPVAR objects or empty doubles');
end
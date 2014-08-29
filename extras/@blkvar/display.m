function F = display(X)
% DISPLAY Overloaded

try
    X = sdpvar(X);
    display(X);
catch
    disp('Incomplete block variable.');
end

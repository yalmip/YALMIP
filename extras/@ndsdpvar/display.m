function F = display(X)
% DISPLAY Overloaded

d = X.dim;
info = num2str(d(1));
for i = 2:length(d)
    info = [info 'x' num2str(d(i))];
end

disp(['Multi-dimensional SDPVAR object ' info])


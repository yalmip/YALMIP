function F = display(X)
% DISPLAY Overloaded

% Author Johan Löfberg
% $Id: display.m,v 1.2 2006-05-17 13:44:06 joloef Exp $

d = X.dim;
info = num2str(d(1));
for i = 2:length(d)
    info = [info 'x' num2str(d(i))];
end

disp(['Multi-dimensional SDPVAR object ' info])


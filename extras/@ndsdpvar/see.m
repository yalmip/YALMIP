function see(X)
% SEE (overloaded)

% Author Johan Löfberg
% $Id: see.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

disp('Constant matrix');disp(' ')
disp(full(getbasematrix(X,0)))
disp('Base matrices');disp(' ')
for i = 1:length(X.lmi_variables);
    disp(full(getbasematrix(X,X.lmi_variables(i))))
    disp(' ')
end;
disp('Used variables');disp(' ')
disp(X.lmi_variables)
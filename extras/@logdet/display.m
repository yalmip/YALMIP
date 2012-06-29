function display(X)
%display           Overloaded

% Author Johan Löfberg 
% $Id: display.m,v 1.2 2007-02-02 09:31:24 joloef Exp $  

P = X.P;
classification = 'Logdet-term ';
[n,m] = size(P);
classification = [classification num2str(n) 'x' num2str(m)];

value = double(X);
if ~isnan(value);
    classification = [classification ' (current value: ' num2str(value) ')'];
end                  
disp(classification);
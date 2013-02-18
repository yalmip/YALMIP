function display(X)
%display           Overloaded

% Author Johan Löfberg 
% $Id: display.m,v 1.3 2007-07-31 13:30:39 joloef Exp $  

n = 0;
for i = 1:length(X.diminOrig)
    n = n + prod(X.diminOrig{i});
end
m = 0;
for i = 1:length(X.dimoutOrig)
    m = m + prod(X.dimoutOrig{i});
end
disp(['Optimizer object with ' num2str(n) ' inputs and ' num2str(m) ' outputs.'])

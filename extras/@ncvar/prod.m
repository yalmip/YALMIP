function product=prod(X)
%PROD (overloaded)

% Author Johan L÷fberg 
% $Id: prod.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   


% n = length(X);
% switch n
%     case 1
%         product = X;
%     case 2
%         x1 = extsubsref(X,1);
%         x2 = extsubsref(X,2);
%         product = x1*x2;
%     case 3
%         x1 = extsubsref(X,1);
%         x2 = extsubsref(X,2);
%         x3 = extsubsref(X,3);
%         product = x1*x2*x3;
%         
%     otherwise
%         m = floor(length(X)/2);
%         x1 = extsubsref(X,1:m);
%         x2 = extsubsref(X,m+1:n);
%         product = prod(x1)*prod(x2);
% end
% 
product = 1;
for i = 1:length(X)
   pick = cell(1,1);pick{1}={i};
   product = product*subsref(X,struct('type','()','subs',pick));
end
function product=prod(X)
%PROD (overloaded)

product = 1;
for i = 1:length(X)
   pick = cell(1,1);pick{1}={i};
   product = product*subsref(X,struct('type','()','subs',pick));
end
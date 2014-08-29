function b = uniquestripped(a)
%UNIQUESTRIPPED  Internal function (version without checkings etc.)

b = sort(a(:)');
b = b(diff([b NaN])~=0);


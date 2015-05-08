function b = uniquestripped(a)
%UNIQUESTRIPPED  Internal function (version without checkings etc.)

b = sort(a(:)');
i = diff([b NaN])~=0;
if ~all(i)
    b = b(i);
end


function y = getbasematrix(x,ind)
%GETBASEMATRIX (overloaded sdpvar/getbasematrix on double)

if ind == 0
    y = x;
else
    y = spalloc(size(x,1),size(x,2),0);
end

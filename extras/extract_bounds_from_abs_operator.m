function LU = extract_bounds_from_abs_operator(LU,extstruct,extvariables,i);

arg = extstruct(i).arg{1};
if is(arg,'lpcone') & isreal(arg)
    vars = getvariables(arg);
    % Absolute value larger than 0
    LU(extvariables(i),1) = max([0 LU(extvariables(i),1)]);
    % Upper bound on linear smaller than abs upper bound
    LU(vars,2) = min([LU(vars,2) LU(extvariables(i),2)]);
    % Lower bound larger than -abs upper
    LU(vars,1) = max([LU(vars,1) -LU(extvariables(i),2)]);
    % Absolute upper bound smaller than max(abs(linear bounds))
    LU(extvariables(i),2) = max(abs(LU(vars,:)));
elseif isreal(arg) & size(arg,1)==1
    absmax = abs(getbase(arg))*[1;max(abs(LU(getvariables(arg),:)),[],2)];
    LU(extvariables(i),1) = max([0 LU(extvariables(i),1)]);
    LU(extvariables(i),2) = min([absmax LU(extvariables(i),2)]);
else
    % At least fix lower bound on absolute
    LU(extvariables(i),1) = max([0 LU(extvariables(i),1)]);
end

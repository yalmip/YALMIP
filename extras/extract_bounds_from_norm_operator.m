function LU = extract_bounds_from_norm_operator(LU,extstruct,extvariables,i);
arg = extstruct(i).arg{1};
if min(size(arg))==1 & length(extstruct(i).arg)>1
    p = extstruct(i).arg{2};
    % Maximize each element conservatively, and use the norm of
    % that worst-case vector
    xmax = abs(getbase(arg))*[1;max(abs(LU(getvariables(arg),:)),[],2)];
    LU(extvariables(i),2) = min([norm(xmax,p) LU(extvariables(i),2)]);
end
% norms are positive...
LU(extvariables(i),1) = max([0 LU(extvariables(i),1)]);


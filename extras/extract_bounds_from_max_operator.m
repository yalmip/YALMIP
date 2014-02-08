function LU = extract_bounds_from_max_operator(LU,extstruct,extvariables,i);
arg = extstruct(i).arg{1};
epi = getvariables(extstruct(i).var);
[M,m] = derivebounds(arg,LU);
LU(epi,1) = max(m);
LU(epi,2) = max(M);


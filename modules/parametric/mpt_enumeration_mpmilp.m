function model = mpt_enumeration_mpmilp(Matrices,options)
% Variable bounds when all binary variables are relaxed
[global_lower,global_upper] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);

Matrices.lb = global_lower;
Matrices.ub = global_upper;
if any(Matrices.lb(end-Matrices.nx+1:end) == Matrices.ub(end-Matrices.nx+1:end))
    model = [];
    return
end

% Enumerate a sufficent set of binary cases
% (exploit SOS and pure binary constraints)
[enums,Matrices] = mpt_enumerate_binary(Matrices);

model = [];
for i = 1:size(enums,2);
    if options.verbose & rem(i,20)==0
        disp(['Binary node ' num2str(i) '/' num2str(size(enums,2))]);
    end
    % Create node problem
    lower = global_lower;
    upper = global_upper;
    lower(Matrices.binary_var_index) = enums(:,i);
    upper(Matrices.binary_var_index) = enums(:,i);
    % Pre-solve, solve and merge
    model = mpt_solvenode(Matrices,lower,upper,Matrices,model,options);
end

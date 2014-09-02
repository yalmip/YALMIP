function sys=double(X)
%DOUBLE Returns current numerical value


solution = yalmip('getsolution');
lmi_variables = X.lmi_variables;
opt_variables = solution.variables;

% Definition of nonlinear variables
mt = yalmip('nonCommutingTable');
mt_local = mt(lmi_variables,2:end);
nterms = sum(mt_local | mt_local,2); 
nonlinears = lmi_variables(nterms > 1);

% Okey, we could not do it really fast...
allextended   = yalmip('extvariables');
if isempty(nonlinears) & isempty(allextended) & all(ismembc(lmi_variables,solution.variables))

    % speed up code for simple linear case
    values = solution.values;
    if isempty(solution.values)
        values = sparse(solution.variables,ones(length(solution.variables),1),solution.optvar,size(mt,1),1);       
        yalmip('setvalues',values);
    else
        if any(isnan(solution.values(lmi_variables(:))))
            values = sparse(solution.variables,ones(length(solution.variables),1),solution.optvar,size(mt,1),1);
            yalmip('setvalues',values);
        end
    end

    sys = X.basis*[1;values(lmi_variables(:))];
    if X.typeflag==10
        sys = eig(full(reshape(sys,X.dim(1),X.dim(2))));
    else
        sys = full(reshape(sys,X.dim(1),X.dim(2)));
    end
    return
end

% All double values
values(size(mt,1),1)=nan;values(:)=nan;
values(solution.variables) = solution.optvar;

% Evaluate the extended operators
if ~isempty(allextended)
    extended_variables = find(ismembc(X.lmi_variables,allextended));
    if ~isempty(extended_variables)
        for i = 1:length(extended_variables)
            extvar = lmi_variables(extended_variables(i));
            extstruct = yalmip('extstruct',extvar);

            for k = 1:length(extstruct.arg)
                if isa(extstruct.arg{k},'sdpvar') | isa(extstruct.arg{k},'constraint')
                    extstruct.arg{k} = double(extstruct.arg{k});
                end
            end

            switch extstruct.fcn
                
                case 'sort'
                    [w,loc] = sort(extstruct.arg{1});
                    if extstruct.arg{2}.isthisloc
                        val = loc(extstruct.arg{2}.i);
                    else
                        val = w(extstruct.arg{2}.i);
                    end
                    
                case 'semivar'
                    % A bit messy, since it not really is a nonlinear
                    % operator. semivar(1) does not return 1, but a new
                    % semivar variable...
                    val = values(getvariables(extstruct.var));

                case 'geomean' % Not 100% MATLAB consistent (Hermitian case differ)
                    val = extstruct.arg{1};
                    if ~any(any(isnan(val)))
                        [n,m] = size(val);
                        if n == m
                            if issymmetric(val)
                                val = max(0,real(det(val)))^(1/n);
                            else
                                val = geomean(val);
                            end
                        else
                            val = geomean(val);
                        end
                    else
                        val = nan;
                    end

                case 'pwf'
                    % Has to be placed here due to the case when
                    % all functions are double, since in this case,
                    % we cannot determine in pwf if we want the double or
                    % create a pw constant function...
                    n = length(extstruct.arg)/2;
                    i = 1;
                    val = nan;
                    while i<=n
                        if min(checkset(extstruct.arg{2*i}))>=0
                            val = extstruct.arg{2*i-1};
                            break
                        end
                        i = i + 1;
                    end

                case 'mpower'      
                    val = extstruct.arg{1};
                case {'max','min'} % Cannot take several inputs, put everything in a vector
                    val = feval(extstruct.fcn,[extstruct.arg{:}]);
                case {'or','and'}
                    try
                        val = feval(extstruct.fcn,extstruct.arg{:});
                    catch
                        val = nan;
                    end

                otherwise 
                    try
                        if ismembc(extvar,allevaluators)
                            % Evaluator variables (used for exp, log, ...)
                            % have an appended argument that we need to
                            % remove. The appended variables is used
                            % internally to convert from general format
                            % f(a'x+c) to f(z), z==a'z+b
                            val = feval(extstruct.fcn,extstruct.arg{1:end-1});
                        else
                            val = feval(extstruct.fcn,extstruct.arg{:});
                        end
                    catch
                        val = nan;
                    end
            end
            values(extvar) = full(val);
        end
    end
end

if ~isempty(nonlinears)
    for i = 1:length(nonlinears)
        index = find(mt(nonlinears(i),2:end));
        index = mt(nonlinears(i),1+index);
        the_product = prod(values(index));
        values(nonlinears(i)) = the_product;
    end
end
sys = X.basis*[1;values(lmi_variables(:))];

if X.typeflag==10
    sys = eig(full(reshape(sys,X.dim(1),X.dim(2))));
else
    sys = full(reshape(sys,X.dim(1),X.dim(2)));
end

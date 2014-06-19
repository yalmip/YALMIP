function xevaled = apply_recursive_evaluation(p,xevaled)

xevaled = xevaled(:)';
for i = 1:length(p.evaluation_scheme)
    switch p.evaluation_scheme{i}.group
        case 'eval'
            xevaled = process_evals(p,xevaled,p.evaluation_scheme{i}.variables);
        case 'monom'
            xevaled = process_monomials(p,xevaled,p.evaluation_scheme{i}.variables);
            xevaled = real(xevaled);
        otherwise
    end
end
xevaled = xevaled(:);

function x = process_monomials(p,x,indicies);
indicies = p.monomials(indicies);
try
    % Do bilinears and quadratics directly
    if max(p.variabletype(indicies))==2
        BilinearIndex = p.variabletype(indicies)==1;
        if any(BilinearIndex)
            Bilinears = indicies(BilinearIndex);
            x(Bilinears) = x(p.BilinearsList(Bilinears,1)).*x(p.BilinearsList(Bilinears,2));
        end
        QuadraticIndex = p.variabletype(indicies)==2;
        if any(QuadraticIndex)
            Quadratics = indicies(QuadraticIndex);
            x(Quadratics) = x(p.QuadraticsList(Quadratics,1)).^2;
        end
    else
        % Mixed stuff. At least do bilinear and quadratics efficiently        
        BilinearIndex = p.variabletype(indicies)==1;
        if any(BilinearIndex)              
            Bilinears = indicies(BilinearIndex);
            x(Bilinears) = x(p.BilinearsList(Bilinears,1)).*x(p.BilinearsList(Bilinears,2));                
            indicies(BilinearIndex) = [];            
        end
        QuadraticIndex = p.variabletype(indicies)==2;
        if any(QuadraticIndex)
            Quadratics = indicies(QuadraticIndex);
            x(Quadratics) = x(p.QuadraticsList(Quadratics,1)).^2;
            indicies(QuadraticIndex)=[];
        end
        
        V = p.monomtable(indicies,:);
        r = find(any(V,1));
        V = V(:,r);
        x(indicies) = prod(repmat(x(r),length(indicies),1).^V,2);                
    end
catch
   
    for i = indicies(:)'
        x(i) =  prod(x.^p.monomtable(i,:),2);
    end
   
end

function x = process_evals(p,x,indicies)
if isfield(p.evalMap{1},'prearg')
    for i = indicies
        arguments = p.evalMap{i}.prearg;
        arguments{1+p.evalMap{i}.argumentIndex} = x(p.evalMap{i}.variableIndex);
        if isequal(arguments{1},'log') & (arguments{1+p.evalMap{i}.argumentIndex}<=0)
            x(p.evalVariables(i)) = -1e4;
        else
            x(p.evalMap{i}.computes(:)) = feval(arguments{:});
        end
    end
else
    for i = indicies
        %arguments = {p.evalMap{i}.fcn,x(p.evalMap{i}.variableIndex)};
        %arguments = {arguments{:},p.evalMap{i}.arg{2:end-1}};
        % Append argument with function name, and remove trailing
        % artificial argument
        arguments =  {p.evalMap{i}.fcn,p.evalMap{i}.arg{1:end-1}};
        arguments{1+p.evalMap{i}.argumentIndex} = x(p.evalMap{i}.variableIndex);
        if isequal(arguments{1},'log') & (arguments{1+p.evalMap{i}.argumentIndex}<=0)
            x(p.evalVariables(i)) = -1e4;  %FIXME DOES NOT WORK
            if length(arguments{2})>1
                disp('Report bug in apply_recursive_evaluation')
            end
        else
            if isfield(p.evalMap{i},'computes')
                x(p.evalMap{i}.computes) = feval(arguments{:});
            else
                x(p.evalVariables(i)) = feval(arguments{:});
            end
        end
    end
end
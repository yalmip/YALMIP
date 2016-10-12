function [model,changed] = convert_sigmonial_to_sdpfun(model)

% Always add this dummy struct
model.high_monom_model = [];

% Assume we don't do anything
changed = 0;

found_and_converted = [];
if any(model.variabletype > 3)
    % Bugger...
    changed = 1;

    % Find a higher order term
    sigmonials = find(model.variabletype == 4);

    model = update_monomial_bounds(model);

    monosig = sigmonials(find(sum(model.monomtable(sigmonials,:) | model.monomtable(sigmonials,:),2)==1));
    if ~isempty(monosig)
        % These are just monomial terms such as x^0.4 etc
        for i = 1:length(monosig)
            variable = find(model.monomtable(monosig(i),:));
            power = model.monomtable(monosig(i),variable);
            model = add_sigmonial_eval(model,monosig(i),variable,power);
            found_and_converted = [found_and_converted;variable power monosig(i)];
        end
    end
end
if any(model.variabletype > 3)
    % Bugger...we have mixed terms such as x/y etc

    % Find a higher order term
    sigmonials = find(model.variabletype == 4);

    for i = 1:length(sigmonials)
        n_old_monoms = size(model.monomtable,1);
        monoms = model.monomtable(sigmonials(i),:);

        % Which variables have fractional or negative powers
        sigs = find((monoms ~= fix(monoms)) | monoms<0);
        powers =  monoms(sigs);
        if ~isempty(found_and_converted)
            % Maybe some of the terms have already been defined as new
            % variables
            for j = 1:length(sigs)
                old_index = findrows(found_and_converted(:,1:2),[sigs(j) powers(j)]);
                if ~isempty(old_index)
                    corresponding_variable = found_and_converted(old_index,3);
                    model.monomtable(sigmonials(i),sigs(j)) = 0;
                    model.monomtable(sigmonials(i),corresponding_variable) = 1;
                    sigs(j)=nan;
                end
            end
        end
        powers(isnan(sigs)) = [];
        sigs(isnan(sigs)) = [];
        if length(sigs) > 0
            % Terms left that haven't been modeled
            model.monomtable(sigmonials(i),sigs) = 0;
            model.monomtable = blkdiag(model.monomtable,speye(length(sigs)));
            model.monomtable(sigmonials(i),n_old_monoms+1:n_old_monoms+length(sigs)) = 1;
            model.variabletype(sigmonials(i)) = 3;
            model.variabletype(end+1:end+length(sigs)) = 0;
            model.c(end+1:end+length(sigs)) = 0;
            model.Q = blkdiag(model.Q,zeros(length(sigs)));
            model.F_struc = [model.F_struc zeros(size(model.F_struc,1),length(sigs))];
            model.lb = [model.lb;-inf(length(sigs),1)];
            model.ub = [model.ub;inf(length(sigs),1)];
            if ~isempty(model.x0)
                model.x0 = [model.x0;model.x0(sigs).^powers(:)];
            end
            for j = 1:length(sigs)                
                model.evalVariables = [model.evalVariables n_old_monoms+j];
                model.isevalVariable(model.evalVariables)=1; 
                if powers(j)==-1
                    model.evalMap{end+1} = inverse_internal2_operator(model,sigs(j),n_old_monoms+j);    
                else
                    model.evalMap{end+1} = power_internal2_operator(model,sigs(j),powers(j));    
                end    
                model.evalMap{end}.properties.domain = [-inf inf];
                model.evalMap{end}.properties.range = [-inf inf];
                model.evalMap{end}.variableIndex = sigs(j);
                model.evalMap{end}.argumentIndex = 1;                
                model.evalMap{end}.computes = n_old_monoms+j;                
                found_and_converted = [found_and_converted;sigs(j) powers(j) n_old_monoms+j];
            end
        end
        if sum(model.monomtable(sigmonials(i),:))<=2
            if nnz(model.monomtable(sigmonials(i),:))==1
                model.variabletype(sigmonials(i)) = 2;
            else
                model.variabletype(sigmonials(i)) = 1;
            end
        end
    end
    
    model = update_eval_bounds(model);
    for i = 1:length(model.evalMap)
        if isequal(model.evalMap{i}.fcn,'power_internal2')
            if isequal(model.evalMap{i}.arg{2},-1)
                if model.lb(model.evalMap{i}.variableIndex) > 0
                    model.evalMap{i}.properties.convexity = 'convex';                    
                    model.evalMap{i}.properties.monotonicity='decreasing';
                    model.evalMap{i}.properties.inverse=@(x)1./x;
                elseif model.ub(model.evalMap{i}.variableIndex) < 0                    
                    model.evalMap{i}.properties.convexity = 'concave';
                    model.evalMap{i}.properties.monotonicity='increasing';
                    model.evalMap{i}.properties.inverse=@(x)1./x;
                end
            end
        end
    end
end

function  model = add_sigmonial_eval(model,monosig,variable,power)
model.evalVariables = [model.evalVariables monosig];
model.isevalVariable(model.evalVariables)=1;  
if power == -1       
    model.evalMap{end+1} = inverse_internal2_operator(model,variable,variable);    
else       
    model.evalMap{end+1} = power_internal2_operator(model,variable,power);
end
model.evalMap{end}.variableIndex = find(model.monomtable(monosig,:));
model.evalMap{end}.argumentIndex = 1;
model.evalMap{end}.computes = monosig;
model.monomtable(monosig,variable) = 0;
model.monomtable(monosig,monosig) = 1;
model.variabletype(monosig) = 0;

% This should not be hidden here....
function [L,U] = power_bound(xL,xU,power)
if xL >= 0
    % This is the easy case
    % we use abs since 0 sometimes actually is -0 but still passes the test
    % above
    if power > 0
        L = abs(xL)^power;
        U = abs(xU)^power;
    else
        L = abs(xU)^power;
        U = abs(xL)^power;
    end
else
    if power < 0 & xU > 0
        % Nasty crossing around zero
        U = inf;
        L = -inf;
    elseif xU < 0
        L = xU^power;
        U = xL^power;
    else
        disp('Not implemented yet')
        error
    end
end

function [L,U] = inverse_bound(xL,xU)
if xL >= 0
    % This is the easy case. We use abs since 0 sometimes actually is -0
    % but still passes the test above
    L = abs(xU)^-1;
    U = abs(xL)^-1;
else
    if xU > 0
        % Nasty crossing around zero
        U = inf;
        L = -inf;
    elseif xU < 0
        L = xU^-1;
        U = xL^-1;
    else
        disp('Not implemented yet')
        error
    end
end

function [Ax, Ay, b] = power_convexhull(xL,xU,power)
fL = xL^power;
fU = xU^power;
dfL = power*xL^(power-1);
dfU = power*xU^(power-1);
if xL<0 & xU>0
    % Nasty crossing
    Ax = [];
    Ay = [];
    b = [];
    return
end

average_derivative = (fU-fL)/(xU-xL);
xM = (average_derivative/power).^(1/(power-1));
if xU < 0
    xM = -xM;
end
fM = xM^power;
dfM = power*xM^(power-1);

if ((power > 1 | power < 0) & (xL >=0)) | ((power < 1 & power > 0) & (xU <=0))
    [Ax,Ay,b] = convexhullConvex(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
else
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
end
if ~isempty(Ax)
    if isinf(Ax(1))
        Ay(1) = 0;
        Ax(1) = -1;
        B(1)  = 0;
    end
end
function [Ax, Ay, b] = inverse_convexhull(xL,xU)
fL = xL^-1;
fU = xU^-1;
dfL = -1*xL^(-2);
dfU = -1*xU^(-2);
if xL<0 & xU>0
    % Nasty crossing
    Ax = [1;-1];
    Ay = [0;0];
    b = [xU;-xL];
    return
end
average_derivative = (fU-fL)/(xU-xL);
xM = (average_derivative/(-1)).^(1/(-1-1));
if xU < 0
    xM = -xM;
end
if ~(xM > xL)
    xM = (xL + xU)/2;
end
fM = xM^(-1);
dfM = (-1)*xM^(-2);

if xL >= 0
    [Ax,Ay,b] = convexhullConvex(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
else
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
end

function df = power_derivative(x,power)
fL = xL^power;
fU = xU^power;
dfL = power*xL^(power-1);
dfU = power*xU^(power-1);
if xL<0 & xU>0
    % Nasty crossing
    Ax = [];
    Ay = [];
    b = [];
    return
end
if power > 1 | power < 0
    [Ax,Ay,b] = convexhullConvex(xL,xU,fL,fU,dfL,dfU);
else
    [Ax,Ay,b] = convexhullConcave(xL,xU,fL,fU,dfL,dfU);
end
if ~isempty(Ax)
    if isinf(Ax(1))
        Ay(1) = 0;
        Ax(1) = -1;
        B(1)  = 0;
    end
end

function  f = inverse_internal2_operator(model,variable,in);
f.fcn = 'inverse_internal2';
f.arg{1} = recover(in);
f.arg{2} = [];
f.properties.bounds = @inverse_bound;
f.properties.convexhull = @inverse_convexhull;
f.properties.derivative = @(x) -1./(x.^2);
f.properties.range = [-inf  inf];
f.properties.domain = [-inf  inf];
flb = 1/model.lb(variable);
fub = 1/model.ub(variable);
if model.lb(variable)>0 | model.ub(variable) < 0
    f.properties.monotonicity = 'decreasing';
    f.properties.inverse = @(x)(1./x);   
    f.properties.range = [min(flb,fub) max(flb,fub)];
end
if model.lb(variable) >= 0
    f.properties.convexity = 'convex';
    f.properties.range = [fub flb];    
elseif model.ub(variable) <= 0
    f.properties.convexity = 'concave';
    f.properties.range = [fub  flb];    
end

function f = power_internal2_operator(model,variable,power);
f.fcn = 'power_internal2';
f.arg{1} = recover(variable);
f.arg{2} = power;
f.arg{3} = [];
f.properties.bounds = @power_bound;
f.properties.convexhull = @power_convexhull;
f.properties.derivative = eval(['@(x) ' num2str(power) '*x.^(' num2str(power) '-1);']);
if even(power)
    f.properties.range = [0 inf];
else
    f.properties.range = [-inf inf];
end
f.properties.domain = [-inf inf];
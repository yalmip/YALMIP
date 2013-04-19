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
                model.evalMap{end+1}.fcn = 'power_internal2';
                model.evalMap{end}.arg{1} = recover(n_old_monoms+j);
                model.evalMap{end}.arg{2} = powers(j);
                model.evalMap{end}.arg{3} = [];
                model.evalMap{end}.variableIndex = sigs(j);
                model.evalMap{end}.argumentIndex = 1;
                model.evalMap{end}.properties.bounds = @power_bound;
                model.evalMap{end}.properties.convexhull = @power_convexhull;
                % TRIAL
                model.evalMap{end}.computes = n_old_monoms+j;
                model.evalMap{end}.properties.derivative = eval(['@(x) ' num2str(powers(j)) '*x.^(' num2str(powers(j)) '-1);']);
                if even(powers(j))
                    model.evalMap{end}.properties.range = [0 inf];
                else
                    model.evalMap{end}.properties.range = [-inf inf];
                end
                model.evalMap{end}.properties.domain = [-inf inf];                                  
                % Save information about this new variable 
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
end

function  model = add_sigmonial_eval(model,monosig,variable,power)
if power == -1
    model.evalVariables = [model.evalVariables monosig];
    model.isevalVariable(model.evalVariables)=1;   
    model.evalMap{end+1}.fcn = 'inverse_internal2';
    model.evalMap{end}.arg{1} = recover(variable);
    model.evalMap{end}.arg{2} = [];
    model.evalMap{end}.variableIndex = find(model.monomtable(monosig,:));
    model.evalMap{end}.argumentIndex = 1;    
    model.evalMap{end}.computes = monosig;
    model.evalMap{end}.properties.bounds = @inverse_bound;
    model.evalMap{end}.properties.convexhull = @inverse_convexhull;
    model.evalMap{end}.properties.derivative = @(x) -1./(x.^2);
    if model.lb(variable)>0 | model.ub(variable) < 0
        model.evalMap{end}.properties.monotonicity = 'decreasing';
        model.evalMap{end}.properties.inverse = @(x)(1./x);
    end
    model.monomtable(monosig,variable) = 0;
    model.monomtable(monosig,monosig) = 1;
    model.variabletype(monosig) = 0;
else
    model.evalVariables = [model.evalVariables monosig];
    model.isevalVariable(model.evalVariables)=1;   
    model.evalMap{end+1}.fcn = 'power_internal2';
    model.evalMap{end}.arg{1} = recover(variable);
    model.evalMap{end}.arg{2} = power;
    model.evalMap{end}.arg{3} = [];
    % Trial
    model.evalMap{end}.computes = monosig;
    model.evalMap{end}.variableIndex = find(model.monomtable(monosig,:));
    model.evalMap{end}.argumentIndex = 1;
    model.evalMap{end}.properties.bounds = @power_bound;
    model.evalMap{end}.properties.convexhull = @power_convexhull;
    if power==fix(power)
        model.evalMap{end}.properties.derivative = eval(['@(x) ' num2str(power) '*x.^(' num2str(power) '-1);']);
    else
        model.evalMap{end}.properties.derivative = eval(['@(x) ' num2str(power) '*max(1e-8,x).^(' num2str(power) '-1);']);
    end
    model.monomtable(monosig,variable) = 0;
    model.monomtable(monosig,monosig) = 1;
    model.variabletype(monosig) = 0;
end

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
    Ax = [];
    Ay = [];
    b = [];
    return
end
average_derivative = (fU-fL)/(xU-xL);
xM = (average_derivative/(-1)).^(1/(-1-1));
if xU < 0
    xM = -xM;
end
fM = xM^(-1);
dfM = (-1)*xM^(-2);

if xL >= 0
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

function [Ax,Ay,b] = fraction_convexhull(xL,xU)
yL = xL(2);
yU = xU(2);
xL = xL(1);
xU = xU(1);

xbar = (xL + xU)/2;
ybar = (yL + yU)/2;

sdpvar x y z

C = [z >= x*(2*(xbar + sqrt(xL*xU)))/((sqrt(xL)+sqrt(xU))^2*ybar) - y*(2*(xbar + sqrt(xL*xU))^2)/((sqrt(xL)+sqrt(xU))^2*ybar^2) + 2*(xbar + sqrt(xL*xU))*sqrt(xL*xU)/(ybar*(sqrt(xL) + sqrt(xU))^2)];
C = [C, z>= xbar/yL - xU/ybar^2*y - (ybar-2*yL)*xU/(ybar*yL)];
C = [C, z>= xbar/yU - xL/ybar^2*y - (ybar-2*yU)*xL/(ybar*yU)];
C = [C, xL<=x<=xU,yL <= y <= yU, z <= xU/yL]
C = [z >= (x*yU - y*xL + xL*yU)/yU^2, z >= (x*yL-y*xU+xU*yL)/yL^2]


for i = 1:100
    for j = 1:100
        xx(i,j) = xL + (xU - xL)*(i-1)/99;
        yy(i,j) = yL + (yU - yL)*(j-1)/99;
        zz(i,j) = xx(i,j)/yy(i,j);
    end
end
Ax = [];
Ay = [];
b = [];

function [L,U] = fraction_bound(xL,xU)
p1 = [xL(1)/xL(2)];
p2 = [xU(1)/xL(2)];
p3 = [xL(1)/xU(2)];
p4 = [xU(1)/xU(2)];
L = min([p1 p2 p3 p4]);
U = max([p1 p2 p3 p4]);

function df =fraction_derivative(xy)
df = [1/xy(2);-xy(1)./xy(2)^2];

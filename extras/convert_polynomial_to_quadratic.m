function [model,changed] = convert_polynomial_to_quadratic(model)

% Assume we don't do anything
changed = 0;

% Are there really any non-quadratic terms?
already_done = 0;
while any(model.variabletype > 2)
    % Bugger...
    changed = 1;

    % Find a higher order term
    polynomials = find(model.variabletype >= 3);
    model = update_monomial_bounds(model,(already_done+1):size(model.monomtable,1));
    already_done = size(model.monomtable,1);
    % Start with the highest order monomial (the bilinearization is not
    % unique, but by starting here, we get a reasonable small
    % bilineared model
    [i,j] = max(sum(abs(model.monomtable(polynomials,:)),2));
    polynomials = polynomials(j);

    powers = model.monomtable(polynomials,:);
    if any(powers < 0)
        model = eliminate_sigmonial(model,polynomials,powers);
    elseif nnz(powers) == 1
        model = bilinearize_recursive(model,polynomials,powers);
    else
        model = bilinearize_recursive(model,polynomials,powers);
    end
end

function   model = eliminate_sigmonial(model,polynomial,powers);
% Silly bug
if isempty(model.F_struc)
    model.F_struc = zeros(0,length(model.c) + 1);
end

% x^(-p1)*w^(p2)
powers_pos = powers;
powers_neg = powers;
powers_pos(powers_pos<0) = 0;
powers_neg(powers_neg>0) = 0;
[model,index_neg] = findoraddlinearmonomial(model,-powers_neg); 
if any(powers_pos)
    [model,index_pos] = findoraddlinearmonomial(model,powers_pos);
else
    index_pos = [];
end
% Now create a new variable y, used to model x^(-p1)*w^(p2)
% We will also add the constraint y*x^p1 = w^p2
model.monomtable(polynomial,:) = 0;
model.monomtable(polynomial,polynomial) = 1;
model.variabletype(polynomial) = 0;
model.high_monom_model = blockthem(model.high_monom_model,[polynomial powers_neg]);;

if ~isempty(model.x0);
    model.x0(polynomial) = prod(model.x0(find(powers_neg))'.^powers_neg);
end
powers = -powers_neg;
powers(polynomial) = 1;
[model,index_xy] = findoraddmonomial(model,powers);
model.lb(index_xy) = 1;
model.ub(index_xy) = 1;

model = convert_polynomial_to_quadratic(model);

function   model = bilinearize_recursive(model,polynomial,powers);
% Silly bug
if isempty(model.F_struc)
    model.F_struc = zeros(0,length(model.c) + 1);
end
% variable^power
if nnz(powers) == 1
    univariate = 1;
    variable = find(powers);
    p1 = floor(powers(variable)/2);
    p2 = ceil(powers(variable)/2);
    powers_1 = powers;powers_1(variable) = p1;
    powers_2 = powers;powers_2(variable) = p2;
else
    univariate = 0;
    variables = find(powers);
    mid = floor(length(variables)/2);
    variables_1 = variables(1:mid);
    variables_2 = variables(mid+1:end);
    powers_1 = powers;
    powers_2 = powers;
    powers_1(variables_2) = 0;
    powers_2(variables_1) = 0;
end

[model,index1] = findoraddlinearmonomial(model,powers_1);
[model,index2] = findoraddlinearmonomial(model,powers_2);
model.monomtable(polynomial,:) = 0;
model.monomtable(polynomial,index1) = model.monomtable(polynomial,index1) + 1;
model.monomtable(polynomial,index2) = model.monomtable(polynomial,index2) + 1;

if index1 == index2
    model.variabletype(polynomial) = 2;
else
    model.variabletype(polynomial) = 1;
end
%model = convert_polynomial_to_quadratic(model);



function [model,index] = findoraddmonomial(model,powers);

if length(powers) < size(model.monomtable,2)
    powers(size(model.monomtable,1)) = 0;
end

index = findrows(model.monomtable,powers);

if isempty(index)
    model.monomtable = [model.monomtable;powers];
    model.monomtable(end,end+1) = 0;
    index = size(model.monomtable,1);
    model.c(end+1) = 0;
    model.Q(end+1,end+1) = 0;
    if size(model.F_struc,1)>0
        model.F_struc(1,end+1) = 0;
    else
        model.F_struc = zeros(0, size(model.F_struc,2)+1);
    end
    bound = powerbound(model.lb,model.ub,powers);
    model.lb(end+1) = bound(1);
    model.ub(end+1) = bound(2);
    if ~isempty(model.x0)
        model.x0(end+1) = 0;
    end
    switch sum(powers)
        case 1
            model.variabletype(end+1) = 0;
        case 2
            if nnz(powers)==1
                model.variabletype(end+1) = 2;
            else
                model.variabletype(end+1) = 1;
            end
        otherwise
            model.variabletype(end+1) = 3;
    end
end

function [model,index] = findoraddlinearmonomial(model,powers);

if sum(powers) == 1
    if length(powers)<size(model.monomtable,2)
        powers(size(model.monomtable,2)) = 0;
    end
    index = findrows(model.monomtable,powers);
    return
end

% We want to see if the monomial x^powers already is modelled by a linear
% variable. If not, we create the linear variable and add the associated
% constraint y == x^powers
index = [];
if ~isempty(model.high_monom_model)
    if length(powers)>size(model.high_monom_model,2)+1
        model.high_monom_model(1,end+1) = 0;
    end
    Z = model.high_monom_model(:,2:end);
    index = findrows(Z,powers);
    if ~isempty(index)
        index = model.high_monom_model(index,1);
        return
    end
end
% OK, we could not find a linear model of this monomial. We thus create a
% linear variable, and add the constraint. Note that it is assumed that the
% nonlinear monomial x^power does exist
[model,index_nonlinear] = findoraddmonomial(model,powers);
model.monomtable(end+1,end+1) = 1;
model.variabletype(end+1) = 0;
model.F_struc = [zeros(1,size(model.F_struc,2));model.F_struc];
model.K.f = model.K.f + 1;
model.F_struc(1,end+1) = 1;
model.F_struc(1,1+index_nonlinear) = -1;
model.c(end+1) = 0;
model.Q(end+1,end+1) = 0;
model.high_monom_model = blockthem(model.high_monom_model,[length(model.c) powers]);;
if ~isempty(model.x0);
    model.x0(end+1) = model.x0(index_nonlinear);
end
model.lb(end+1) = model.lb(index_nonlinear);
model.ub(end+1) = model.ub(index_nonlinear);
index = length(model.c);

function z = initial(x0,powers)
z = 1;
vars  = find(powers);
for k = 1:length(vars)
    z = z * x0(vars(k))^powers(vars(k));
end

function A = blockthem(A,B)
n = size(A,2);
m = size(B,2);
A = [A zeros(size(A,1),max([0 m-n]));B zeros(size(B,1),max([0 n-m]))];
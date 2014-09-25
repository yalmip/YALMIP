function F = projection(F,x,method,alreadyprojected)
% projection  Projects polytopic set object (Requires the Multi-parametric Toolbox).
%
% Fproj     = projection(F,x)
%
% F      : Polytopic set object
% x      : Variables to project on
% method : See HELP PROJECTION

if nargin<2
    error('Not enough input arguments.')
end

if nargin < 4    
    [F,dummy] = expandmodel(F,[]);
elseif ~alreadyprojected
    [F,dummy] = expandmodel(F,[]);
end

F = flatten(F);

f = [];
e = [];
for i = 1:length(F)
    fi =  F.clauses{i}.data;
    if  F.clauses{i}.type==2        
        f = [f;fi(:)];
    else
        e = [e;fi(:)];       
    end
end

if ~islinear(F)
    error('Only linear element-wise inequalities can be projected')
end

allConstraints = [e;f];
B = full(getbase(allConstraints));

P = polytope(-B(:,2:end),B(:,1));
y_vars = getvariables(allConstraints);
x_vars = getvariables(x);

x_vars = [];
for i = 1:length(x)
    x_vars = [x_vars getvariables(x(i))];
end

if ~isempty(setdiff(x_vars,y_vars))
    error('Huh? Project on what! Some of those variables are not part of the original polyhedron');
end

x_in_y = find(ismember(y_vars,x_vars));

% Okay, first length(e) rows ocrresponds to equality constraint
if length(e)>0
    m = length(e);
    beq = B(1:m,1);
    bin = B(m+1:end,1);
    A = -B(:,2:end);
    keepvar = find(ismember(y_vars,x_vars));
    zvar    = find(~ismember(y_vars,x_vars));
    Aeqx = A(1:m,keepvar);
    Aeqz = A(1:m,zvar);
    Ainx = A(m+1:end,keepvar);
    Ainz = A(m+1:end,zvar);
    P = polytope([Ainx-Ainz*pinv(Aeqz)*Aeqx -Ainz*null(Aeqz)],bin-Ainz*pinv(Aeqz)*beq);
    if isfulldim(P)
        P = projection(P,1:length(x));
        H = get(P,'H');
        K = get(P,'K');

        [ii,jj] = sort(x_vars);
        x = recover(x_vars(jj));
        F = (H(:,jj)*x <= K);
    else
        F=(sum(x) <= -1) + (sum(x) >= 1);
    end
else

    if nargin == 2
        P = projection(P,x_in_y);
    elseif ~isempty(method)
        Options.projection=method;
        P = projection(P,x_in_y,Options);
    else
        P = projection(P,x_in_y);
    end

    H = get(P,'H');
    K = get(P,'K');

    [ii,jj] = sort(x_vars);
    x = recover(x_vars(jj));
    F = (H(:,jj)*x <= K);
end

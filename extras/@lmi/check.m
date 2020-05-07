function [pres,dres] = check(F)
% CHECK(F)  Displays/calculates constraint residuals on constraint F
%
% [pres,dres] = check(F)
%
% pres : Primal constraint residuals
% dres : Dual constraint residuals
%
% If no output argument is supplied, tabulated results are displayed
%
% Primal constraint residuals are calculated as:
%
%  Semidefinite constraint F(x)>=0 : min(eig(F))
%  Element-wise constraint F(x)>=0 : min(min(F))
%  Equality constraint F==0        : -max(max(abs(F)))
%  Second order cone t>=||x||      : t-||x||
%  Exponential cone constraint     : x(3) - exp(x(1)/x(2))*x(2)
%  Integrality constraint on x     : max(abs(x-round(x)))
%  Sum-of-square constraint        : Minus value of largest (absolute value) coefficient 
%                                   in the polynomial p-v'*v
%
% Dual constraints are evaluated similarily.
%
%    See also   SOLVESDP, SOLVESOS, SOSD, DUAL

% Check if solution avaliable
currsol = evalin('caller','yalmip(''getsolution'')');
if isempty(currsol)
    disp('No solution available.')
    return
end

F = flatten(F);
nlmi = length(F.LMIid);
spaces = ['                                    '];
if (nlmi == 0)
    if nargout == 0
        disp('empty LMI')
    else
        pres = [];
        dres = [];
    end 
    return
end

lmiinfo{1} = 'Matrix inequality';
lmiinfo{2} = 'Elementwise inequality';
lmiinfo{3} = 'Equality constraint';
lmiinfo{4} = 'Second order cone constraint';
lmiinfo{5} = 'Rotated Lorentz constraint';
lmiinfo{7} = 'Integer constraint';
lmiinfo{8} = 'Binary constraint';
lmiinfo{9} = 'KYP constraint';
lmiinfo{10} = 'Eigenvalue constraint';
lmiinfo{11} = 'SOS constraint';
lmiinfo{15} = 'Uncertainty declaration';
lmiinfo{54} = 'Vectorized second order cone constraints';
lmiinfo{55} = 'Complementarity constraint';
lmiinfo{56} = 'Logic constraint';
header = {'ID','Constraint','Primal residual','Dual residual','Tag'};

if nargout==0
    disp(' ');
end

% Checkset is very slow for multiple SOS
% The reason is that REPLACE has to be called
% for every SOS. Instead, precalc on one vector
p=[];
p_dim = [];
ParVar=[];
soscount=1;
for j = 1:nlmi
    if F.clauses{j}.type==11
        pi = F.clauses{j}.data;   
        [v,ParVari] = sosd(pi);
        if isempty(v)
            p=[p;0];
            p_dim = [p_dim size(pi,1)];
        else
            p=[p;reshape(pi,[],1)];
            p_dim = [p_dim size(pi,1)];
            ParVar=unique([ParVar(:);ParVari(:)]);
        end
    end
end
if ~isempty(ParVar)
    ParVar = recover(ParVar);
    pp = replace(p,ParVar,double(ParVar));
    top = 1;
    clear p
    for j = 1:length(p_dim)
        p{j} = reshape(pp(top:top+p_dim(j)^2-1),p_dim(j),p_dim(j));
        top = top + p_dim(j)^2;
    end
end 

for j = 1:nlmi
    constraint_type = F.clauses{j}.type;
    if constraint_type~=11 && constraint_type~=56
        F0 = double(F.clauses{j}.data);
    end
    if ~((constraint_type == 56) || (constraint_type==11)) && any(isnan(F0(:)))
        res(j,1) = NaN;
    else
        switch F.clauses{j}.type
        case {1,9}
            if isa(F0,'intval')
                res(j,1) = full(min(inf_(veigsym(F0))));
            else
                res(j,1) = full(min(real(eig(F0))));                
            end
        case 2
            if isa(F0,'intval')
                res(j,1) = full(min(min(inf_(F0))));
            else
                res(j,1) = full(min(min(F0)));
            end
        case 3
            res(j,1) = -full(max(max(abs(F0))));
        case 4
            res(j,1) = full(F0(1)-norm(F0(2:end)));
        case 5
            res(j,1) = full(2*F0(1)*F0(2)-norm(F0(3:end))^2);
        case 7
            res(j,1) = -full(max(max(abs(F0-round(F0)))));            
        case 8
            res(j,1) = -full(max(max(abs(F0-round(F0)))));
            res(j,1) = min(res(j,1),-(any(F0>1) | any(F0<0)));
        case 21
            res(j,1) = F0(3) - F0(2)*exp(F0(1)/F0(2));
        case 54
            res(j,1) = inf;
            for k = 1:size(F0,2)
                res(j,1) = min(res(j,1),full(F0(1,k)-norm(F0(2:end,k))));
            end                
        case 11
                 
            [v,ParVar] = sosd(F.clauses{j}.data);
            if isempty(v)
                res(j,1)=nan;
            else
                h = p{soscount}-v'*v;
                res(j,1) = full(max(max(abs(getbase(h)))));
            end
            soscount=soscount+1;
            
            case 56
                 res(j,1) = logicSatisfaction(F.clauses{j}.data);
            
        otherwise
            res(j,1) = nan;
        end
    end
    
    % Get the internal index   
    LMIid = F.LMIid(j);
    dual  = yalmip('dual',LMIid);
    if isempty(dual) | any(isnan(dual))
        resdual(j,1) = NaN;
    else
        switch F.clauses{j}.type
            case {1,9}
                resdual(j,1) = min(eig(dual));
            case 2
                resdual(j,1) = min(min(dual));
            case 3
                resdual(j,1) = 0;
            case 4
                resdual(j,1) = dual(1)-norm(dual(2:end));
            case 5
                resdual(j,1) = 2*dual(1)*dual(2)-norm(dual(3:end))^2;
            case 7
                resdual(j,1) = nan;
            case 54
                resdual(j,1) = inf;
                for k = 1:size(dual,2)
                    resdual(j,1) = min(resdual(j,1),full(dual(1,k)-norm(dual(2:end,k))));
                end
            otherwise
                gap = nan;
        end
    end
    
    if nargout==0
        data{j,1} = ['#' num2str(j)];        
        data{j,2} = lmiinfo{F.clauses{j}.type};
        data{j,3} = res(j,1);
        data{j,4} = resdual(j,1);
        data{j,5} = F.clauses{j}.handle;
        
        
        if ~islinear(F.clauses{j}.data)
            if is(F.clauses{j}.data,'sigmonial')
                classification = 'sigmonial';
            elseif is(F.clauses{j}.data,'bilinear')
                classification = 'bilinear';
            elseif is(F.clauses{j}.data,'quadratic')
                classification = 'quadratic';
            else
                classification = 'polynomial';
            end
            data{j,2} = [data{j,2} ' (' classification ')'];
        end
end
end

if nargout>0
    pres = res;
    dres = resdual;
else
    keep = ones(1,5);
    
    if length([data{:,5}])==0
        keep(5) = 0;
    end
    
    header = {header{:,find(keep)}};
    temp = {data{:,find(keep)}};
    data = reshape(temp,length(temp)/nnz(keep),nnz(keep));

    post{1} = 'A primal-dual optimal solution would show non-negative residuals.';
    post{end+1} = 'In practice, many solvers converge to slightly infeasible';
    post{end+1} = 'solutions, which may cause some residuals to be negative.';
    post{end+1} = 'It is up to the user to judge the importance and impact of';
    post{end+1} = 'slightly negative residuals (i.e. infeasibilities)';
    post{end+1} = 'https://yalmip.github.io/command/check/';
    post{end+1} = 'https://yalmip.github.io/faq/solutionviolated/';
    yalmiptable('',header,data,[],post);
    disp(' ');
end

function res = logicSatisfaction(clause);

if isequal(clause{1},'ismember')
    if isa(clause{2},'sdpvar') &  isa(clause{3},'double')
        res = -min(sum(abs(value(clause{2}) - clause{3}),1))
        return
    end
end
a = clause{2};
b = clause{3};
if isa(a,'sdpvar')
    aval = double(a);
    if is(a,'binary')
        atruth = aval == 1;
    else
        atruth = aval>=0;
    end
elseif isa(a,'lmi') | isa(a,'constraint')
    aval = check(lmi(a));
    atruth = aval >= 0;
end
if isa(b,'sdpvar')
    bval = double(b);
    if is(b,'binary')
        btruth = bval == 1;
    else
        btruth = bval>=0;
    end
elseif isa(b,'lmi') | isa(b,'constraint')
    bval = check(lmi(b));
    btruth = bval >= 0;
end
if isnan(aval) | isnan(bval)
    res = nan;
else
    switch clause{1}
        case 'implies'
            if all(btruth >= atruth)
                res = 1;
            else
                res = -1;
            end
        case 'iff'
            if all(btruth == atruth);
                res = 1;
            else
                res = -1;
            end
        otherwise
            res = nan;
    end
end
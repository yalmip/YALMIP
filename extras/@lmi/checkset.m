function [pres,dres] = checkset(F)
% CHECKSET(F)  Displays/calculates constraint residuals on constraint F
%
% [pres,dres] = checkset(F)
%
% pres : Primal constraint residuals
% dres : Dual constraint residuals
%
% If no output argument is supplied, tabulated results are displayed
%
% Primal constraint residuals are calculated as:
%
%  Semidefinite constraint F(x)>0 : min(eig(F))
%  Element-wise constraint F(x)>0 : min(min(F))
%  Equality constraint F==0       : -max(max(abs(F)))
%  Second order cone t>||x||      : t-||x||
%  Integrality constraint on x    : max(abs(x-round(x)))
%  Rank constraint rank(X) < r     : r-rank(X)
%  Sum-of-square constraint       : Minus value of largest (absolute value) coefficient 
%                                   in the polynomial p-v'*v
%
% Dual constraints are evaluated similarily.
%
%    See also   SET, SOLVESDP, SOLVESOS, SOSD, DUAL

% Author Johan Löfberg
% $Id: checkset.m,v 1.19 2008-02-21 15:32:09 joloef Exp $

% Check if solution avaliable
currsol = evalin('caller','yalmip(''getsolution'')');
if isempty(currsol)
    disp('No solution available.')
    return
end

nlmi = size(F.clauses,2);
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
header = {'ID','Constraint','Type','Primal residual','Dual residual','Tag'};

if nargout==0
    disp(' ');
end

% Checkset is very slow for multiple SOS
% The reason is that REPLACE has to be called
% for every SOS. Instead, precalc on one vector
p=[];
ParVar=[];
soscount=1;
for j = 1:nlmi
    if F.clauses{j}.type==11
        pi = F.clauses{j}.data;   
        [v,ParVari] = sosd(pi);
        if isempty(v)
            p=[p;0];
        else
            p=[p;pi];
            ParVar=unique([ParVar(:);ParVari(:)]);
        end
    end
end
if ~isempty(ParVar)
    ParVar = recover(ParVar);
    p = replace(p,ParVar,double(ParVar));
end 

for j = 1:nlmi
    constraint_type = F.clauses{j}.type;
    if constraint_type~=11
        F0 = double(F.clauses{j}.data);
    end
    if (constraint_type~=11) & any(isnan(F0(:)))
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
        case 54
            res(j,1) = inf;
            for k = 1:size(F0,2)
                res(j,1) = min(res(j,1),full(F0(1,k)-norm(F0(2:end,k))));
            end                
        case 11
            if 0
                p = F.clauses{j}.data;          
                [v,ParVar] = sosd(p);
                if ~isempty(ParVar)
                    ParVar = recover(ParVar);
                    p = replace(p,ParVar,double(ParVar));
                end
                if isempty(v)
                    res(j,1)=nan;
                else
                    h = p-v'*v;
                    res(j,1) = full(max(max(abs(getbase(h)))));
                end
            else
                %p = F.clauses{j}.data;          
                [v,ParVar] = sosd(F.clauses{j}.data);
                if isempty(v)
                    res(j,1)=nan;
                else
                    h = p(soscount)-v'*v;
                    res(j,1) = full(max(max(abs(getbase(h)))));
                end
                soscount=soscount+1;
            end
            
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
            resdual(j,1) = -max(max(abs(dual)));
        case 4
            resdual(j,1) = dual(1)-norm(dual(2:end));
        case 5
            resdual(j,1) = 2*dual(1)*dual(2)-norm(dual(3:end))^2;
        case 7
            resdual(j,1) = nan;
        otherwise
            gap = nan;
        end
    end
%     if isempty(dual) | any(isnan(dual)) |  any(isnan(F0))
%         gap = NaN;
%     else
%         switch F.clauses{j}.type
%         case {1,9}
%             gap = trace(F0*dual);
%         case {2,3}
%             gap = F0(:)'*dual(:);
%         case 4
%             gap = nan;
%         case 5
%             gap = nan;
%         case 7
%             gap = nan;
%         otherwise
%             gap = nan;
%         end
%     end
    
    if nargout==0
        data{j,1} = ['#' num2str(j)];
        data{j,2} = F.clauses{j}.symbolic;
        data{j,3} = lmiinfo{F.clauses{j}.type};
        data{j,4} = res(j,1);
        data{j,5} = resdual(j,1);
%        data{j,6} = gap;
        data{j,6} = F.clauses{j}.handle;
        
        
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
            data{j,3} = [data{j,3} ' (' classification ')'];
        end
end
end

if nargout>0
    pres = res;
    dres = resdual;
else
    if length([data{:,6}])==0
        header = {header{:,1:5}};
        temp = {data{:,1:5}};
        data = reshape(temp,length(temp)/5,5);
    end
    
    if all(isnan(resdual))
        header = {header{:,[1 2 3 4]}};
        temp = {data{:,[1 2 3 4]}};
        data = reshape(temp,length(temp)/4,4);
    end
    yalmiptable('',header,data)
    disp(' ');
end
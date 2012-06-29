function [lb,ub,redundant,psstruct,infeasible] = milppresolve(A,b,lb,ub,integer_variables,binary_variables,changed_bounds);
%MILPPRESOLVE Internal function for presolving MILPs

% Author Johan Löfberg
% $Id: milppresolve.m,v 1.3 2005-10-05 14:54:55 joloef Exp $

% Simple bound pre-processing (paper by Savelsbergh)
% No code optimization at all

if nargin < 5
    integer_variables = [];
end
if nargin < 6
    binary_variables = [];
end
if nargin < 7
    changed_bounds = [];
end

ubnew = ub;
lbnew = lb;
goon = 1;

AL0A  = (A>0).*A;
AG0A  = (A<0).*A;

At = A';
use_indicies=ones(length(b),1);
used = full(any(A(:,find(changed_bounds)),2));
isbinary  = ismembc(1:length(lb),binary_variables);
isinteger = ismembc(1:length(lb),binary_variables);

goon = all(lb<=ub);
infeasible = ~goon;
if ~infeasible
    
    bi_up = AL0A*ub+AG0A*lb;
    bi_dn = AL0A*lb+AG0A*ub;

    while goon

        bminusbdn=b-bi_dn;
        for i = find(use_indicies & used)'

            [ii,jj,kk] = find(At(:,i));

            Cp = ii(kk>0);
            if ~isempty(Cp)
                new1=lbnew(Cp)+bminusbdn(i)./At(Cp,i);
                ubnew(Cp)=min(ubnew(Cp),new1);
            end

            Cm = ii(kk<0);
            if ~isempty(Cm)
                new2=ubnew(Cm)+bminusbdn(i)./At(Cm,i);
                lbnew(Cm)=max(lbnew(Cm),new2);
            end

            if any(lbnew>ubnew)
                infeasible = 1;
                break
            end

        end

        if ~infeasible

            lbnew(integer_variables) = round(lbnew(integer_variables)+0.49999);
            ubnew(integer_variables) = round(ubnew(integer_variables)-0.49999);

            lbnew(binary_variables) = round(lbnew(binary_variables)+0.49999);
            ubnew(binary_variables) = round(ubnew(binary_variables)-0.49999);

            goon = (~all((lb==lbnew) & (ub==ubnew))) & all(lbnew<=ubnew);

            used = (lb~=lbnew) | (ub~=ubnew);
            used = full(any(A(:,find(used)),2));

            lb = lbnew;
            ub = ubnew;                                  
            
            bi_up = AL0A*ub+AG0A*lb;
            bi_dn = AL0A*lb+AG0A*ub;            
            redundant = find(bi_up<=b);   

            use_indicies = use_indicies & bi_up>b;
        else
            goon = 0;
        end
    end
    redundant = find(bi_up<=b);
else
    redundant=[];
end
psstruct.AL0A = AL0A;
psstruct.AG0A = AG0A;



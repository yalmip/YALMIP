function [lb,ub,redundant,psstruct,infeasible] = tightenbounds(A,b,lb,ub,integer_variables,binary_variables,changed_bounds);
%TIGHTENBOUNDS Internal function to perform bound tightening

% Author Johan Löfberg
% $Id: tightenbounds.m,v 1.1 2006-03-30 13:56:54 joloef Exp $

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
Cpp = cell(1,size(At,2));
Cmm = cell(1,size(At,2));
r1 = cell(1,size(At,2));
for i = 1:size(At,2)
    [ii,jj,kk] = find(At(:,i));
    Cpp{i} = ii(kk>0);
    Cmm{i} = ii(kk<0);
end

iter1 = 1;
if ~infeasible & (iter1 < 10)
    iter2 = 1;
    while goon & iter2 < 10

        bi_up = AL0A*ub+AG0A*lb;
        bi_dn = AL0A*lb+AG0A*ub;
% 
         if ~isempty(find(bi_dn>b))
             infeasible = 1;
             break;
         end
%         
        bminusbdn=b-bi_dn;
        for i = find(use_indicies & used)'


            Cp = Cpp{i};
            if ~isempty(Cp)
                r = At(Cp,i);              
                new1=lbnew(Cp)+bminusbdn(i)./r;
                ubnew(Cp) = min(ubnew(Cp),new1);
            end
            
            Cm = Cmm{i};
            if ~isempty(Cm)
                r = At(Cm,i);
                new2=ubnew(Cm)+bminusbdn(i)./r;
                lbnew(Cm) = max(lbnew(Cm),new2);
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

            use_indicies = use_indicies & bi_up>b;
        else
            goon = 0;
        end
        iter2 = iter2 + 1;
    end
    iter1 = iter1 + 1;
    redundant = find(bi_up<=b-0.001);
else
    redundant=[];
end
psstruct.AL0A = AL0A;
psstruct.AG0A = AG0A;



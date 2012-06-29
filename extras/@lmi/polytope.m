function [P,x] = polytope(X,options)
% polytope  Converts set object to polytope object        
%
% P     = polytope(F)
% [P,x] = polytope(F)
%
% P : polytope object (Requires the Multi-parametric Toolbox)
% x : sdpvar object defining the variables in the polytope P.H*x<P.K
% F : set-object with linear inequalities

% Author Johan Löfberg
% $Id: polytope.m,v 1.3 2005-02-04 10:10:27 johanl Exp $


%[model,recoverdata,diagnostic,p] = export(X,[],[],[],[],0);
if nargin < 2
    options = sdpsettings;
elseif isempty(options)
    options = sdpsettings;
end
    
[p,recoverdata,solver,diagnostic,F] = compileinterfacedata(X,[],[],[],options,0);

if p.K.q(1) > 0 | p.K.s(1) > 0 | any(p.variabletype)
  error('Polytope can only be applied to MILP-representable constraints.')
end

if isempty(p.binary_variables) & isempty(p.integer_variables)
    P = polytope(-p.F_struc(:,2:end),p.F_struc(:,1));
    x = recover(p.used_variables);
else
    
    nBin = length(p.binary_variables);
    [pBinary,removeEQ,removeLP] = extractOnly(p,p.binary_variables);
    p.F_struc = [p.F_struc(removeEQ,:);p.F_struc(p.K.f+removeLP,:)];
    p.K.f = length(removeEQ);
    p.K.l = length(removeLP);
    
    if p.K.f > 0
        disp('MPT does not support polytopes with empty interior')
        disp('Note that these equality constraints might have been generated internally by YALMIP')
        error('Functionality not yet supported')        
    end
    
    P = [];
    for i = 0:2^nBin-1;
        comb = dec2decbin(i,nBin);
        if checkfeasiblefast(pBinary,comb(:),1e-6)
            pi = p;
            H = -p.F_struc(:,2:end);% Hx < K
            K = p.F_struc(:,1);
            K = K-H(:,p.binary_variables)*comb(:);
            H(:,p.binary_variables)=[];
            P = [P polytope(H,K)];
        end
    end
end
%     f = [];
%     X = expandmodel(X,[]);
%     for i = 1:length(X)
%         if  X.clauses{i}.type==2
%             fi =  X.clauses{i}.data;
%             f = [f;fi(:)];
%         end
%     end
%     B = full(getbase(f));
%     p = polytope(-B(:,2:end),B(:,1));




function pLP = extractLP(p);
pLP = p;
pLP.F_struc = pLP.F_struc(1:p.K.f+p.K.l,:);
pLP.K.q = 0;
pLP.K.s = 0;

function [pRed,removeEQ,removeLP] = extractOnly(p,these);
pRed = p;
p.F_struc(:,1+these) = 0;

removeEQ = find(any(p.F_struc(1:pRed.K.f,2:end),2));
removeLP = find(any(p.F_struc(1+pRed.K.f:end,2:end),2));
pRed.F_struc(pRed.K.f+removeLP,:)=[];
pRed.F_struc(removeEQ,:)=[];
pRed.K.f = pRed.K.f - length(removeEQ);
pRed.K.l = pRed.K.l - length(removeLP);
pRed.F_struc = pRed.F_struc(:,[1 1+these]);
pRed.lb = pRed.lb(these);
pRed.ub = pRed.ub(these);
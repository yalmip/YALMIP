function [Fconv,info] = nonlineartocone(F)
%NONLINEARTOCONE Convert nonlinear constraint to second order cone

% Author Johan Löfberg 
% $Id: nonlineartocone.m,v 1.1 2006-08-10 18:00:21 joloef Exp $   

if islinear(F)
    Fconv = F;
    info = 0;
    return
end

[n,m]=size(F);
if (n*m==1)
    Fconv = [];
    % Involved in polynomial constraints
    sqrList = yalmip('nonlinearvariables');
    h = sort(unique(sqrList));
    % % SDP relaxation
    for i = 1:size(sqrList,1)
        left_index = find(h==sqrList(i,1));
        x1_index = find(h==sqrList(i,2));
        x2_index = find(h==sqrList(i,3));
        Z{i}=zeros(length(h),length(h));
        if x1_index==x2_index
            Z{i}(x2_index,x1_index) = 1;
        else
            Z{i}(x1_index,x2_index) = 0.5;
            Z{i}(x2_index,x1_index) = 0.5;
        end
    end
    vars = getvariables(F);
    base = getbase(F);
    nonlinvars = find(ismember(vars,sqrList(:,1)));
    linvars = find(~ismember(vars,sqrList(:,1)));
    % Construct quadratic 
    Q = zeros(length(h));
    for i = 1:length(nonlinvars)
        indexinlist = find(vars(nonlinvars(i))==sqrList(:,1));
        indexinh = find(vars(nonlinvars(i))==h);
        Q = Q-Z{indexinlist}*base(1+find(h(indexinh)==vars));
    end
    used = find(any(Q));
    [B,r] = chol(Q(used,used));
    if r==0
        linear = base(1);
        if ~isempty(linvars)
            linear = linear+base(1+find(ismember(linvars,vars)))*recover(linvars);
        end
        ConesAxb=[2*B*recover(h(used));1-linear];
        Conescxd=1+linear;
        Fconv = cone(ConesAxb,Conescxd);
        info = 0;
    else
        Fconv = F;
        if nargout == 2
            info = 1;
        else
            warning('Cannot re-write to second order cone constraint');
        end
    end
else
    if nargout == 2
        Fconv = F;
        info = 1;
    else
        error('nonlineartocone can only be applied to scalar inequalities')
    end
end

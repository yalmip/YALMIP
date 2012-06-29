function [g,geq,dg,dgeq] = sqplab_con(x)

global SQPLABDATA

% Early bail for linear problems
if SQPLABDATA.linearconstraints
    if ~isempty(SQPLABDATA.A)
        g = SQPLABDATA.A*x;
        dg = SQPLABDATA.A';
    else
        g = [];
        dg = [];
    end
    if ~isempty(SQPLABDATA.Aeq)
        geq = SQPLABDATA.Aeq*x - SQPLABDATA.beq;
        dgeq = SQPLABDATA.Aeq';
    else
        geq = [];
        dgeq = [];
    end
    return
end

xevaled = zeros(1,length(SQPLABDATA.interfacedata.c));
xevaled(SQPLABDATA.linearindicies) = x;

% Experimental support for arbitrary functions
% nonlinear expressions inside sin() exp() etc
if SQPLABDATA.evalinconstraint
    
    if ~isempty(SQPLABDATA.bilinears)
        xevaled(SQPLABDATA.bilinears(:,1)) = xevaled(SQPLABDATA.bilinears(:,2)).*xevaled(SQPLABDATA.bilinears(:,3));
    else
        pp = SQPLABDATA.monomtable(SQPLABDATA.nonlinearindicies,:);       
        xevaled(SQPLABDATA.nonlinearindicies) = prod(repmat(xevaled,length(SQPLABDATA.nonlinearindicies),1).^SQPLABDATA.monomtable(SQPLABDATA.nonlinearindicies,:),2);
    end

    for i = 1:length(SQPLABDATA.interfacedata.evalMap)
        arguments = SQPLABDATA.evalMap{i}.prearg;
        arguments{2} = xevaled(SQPLABDATA.interfacedata.evalMap{i}.variableIndex);
        xevaled(SQPLABDATA.interfacedata.evalVariables(i)) = feval(arguments{:});
    end
end

if ~isempty(SQPLABDATA.bilinears)
    xevaled(SQPLABDATA.bilinears(:,1)) = xevaled(SQPLABDATA.bilinears(:,2)).*xevaled(SQPLABDATA.bilinears(:,3));
else
    xevaled(SQPLABDATA.nonlinearindicies) = prod(repmat(xevaled,length(SQPLABDATA.nonlinearindicies),1).^SQPLABDATA.monomtable(SQPLABDATA.nonlinearindicies,:),2);
end

if SQPLABDATA.nonlinearinequalities
    g = SQPLABDATA.interfacedata.Anonlinineq*xevaled(:)-0*SQPLABDATA.interfacedata.bnonlinineq;
else
    g = [];
end

if SQPLABDATA.nonlinearequalities
    geq = SQPLABDATA.interfacedata.Anonlineq*xevaled(:)-SQPLABDATA.interfacedata.bnonlineq;
else
    geq = [];
end

K = SQPLABDATA.interfacedata.K;
top = 1;
if K.q(1) > 0
    for i = 1:length(K.q)
        Axcd = SQPLABDATA.interfacedata.F_struc(top:top+K.q(i)-1,:)*[1;xevaled(:)];
        g = [g;-(Axcd(1)^2-norm(Axcd(2:end),2)^2)];
        top = top + K.q(i);
    end
end
if K.s(1) > 0
    for i = 1:length(K.s)
        CminusA = SQPLABDATA.interfacedata.F_struc(top:top+K.s(i)^2-1,:)*[1;xevaled(:)];
        CminusA = reshape(CminusA,K.s(i),K.s(i));
        [R,p] = chol(CminusA);
        if p
            g = [g;-min(eig(CminusA))];
        else
            g = [g;-log(det(CminusA))];
        end
        top = top + K.s(i)^2;
    end
end

dg = [];
dgeq = [];
if SQPLABDATA.SimpleNonlinearConstraints
    dg = [];
    allA = [SQPLABDATA.interfacedata.Anonlineq;SQPLABDATA.interfacedata.Anonlinineq];
    for i = 1:length(SQPLABDATA.linearindicies)
        xevaled = zeros(1,length(SQPLABDATA.interfacedata.c));
        xevaled(SQPLABDATA.linearindicies) = x;
        mt = SQPLABDATA.monomtable;
        oldpower = mt(:,SQPLABDATA.linearindicies(i));
        mt(:,SQPLABDATA.linearindicies(i)) = mt(:,SQPLABDATA.linearindicies(i))-1;
        xevaled = prod(repmat(xevaled,size(mt,1),1).^mt,2);
        xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
        dg = [dg allA*xevaled'];
    end
    dgeq = dg(1:size(SQPLABDATA.interfacedata.Anonlineq,1),:)';
    dg = dg(size(SQPLABDATA.interfacedata.Anonlineq,1)+1:end,:)';
end

%
% dg = [];
% dgeq = [];
% if SQPLABDATA.SimpleNonlinearConstraints
%     dg = [];
%     allA = [SQPLABDATA.interfacedata.Anonlineq;SQPLABDATA.interfacedata.Anonlinineq];
%      mt = SQPLABDATA.monomtable;
%     xe = zeros(1,length(SQPLABDATA.interfacedata.c));
%     xe(SQPLABDATA.linearindicies) = x;
%     xx = xe;
%     xe = repmat(xe,size(mt,1),1).^mt;
%
%     for i = 1:length(SQPLABDATA.linearindicies)
%        % xevaled = zeros(1,length(SQPLABDATA.interfacedata.c));
%        % xevaled(SQPLABDATA.linearindicies) = x;
%
%         oldpower = mt(:,SQPLABDATA.linearindicies(i));
%         newpower = oldpower-1;
% %        mt(:,SQPLABDATA.linearindicies(i)) = mt(:,SQPLABDATA.linearindicies(i))-1;
%         xevaled = xe;xevaled(SQPLABDATA.linearindicies(i),:) = xx(SQPLABDATA.linearindicies(i),:)
%         xevaledSQPLABDATA.linearindicies(i)) = xevaled(SQPLABDATA.linearindicies(i),:).^newpower';
%         xevaled = prod(xevaled,2);
%         xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
%         dg = [dg allA*xevaled'];
%     end
%     dgeq = dg(1:size(SQPLABDATA.interfacedata.Anonlineq,1),:)';
%     dg = dg(size(SQPLABDATA.interfacedata.Anonlineq,1)+1:end,:)';
%     full(dgeq)
% end
%




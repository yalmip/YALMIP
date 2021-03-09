function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_primal(model);

prob.c = model.c;
if ~isempty(model.F_struc)
    prob.a = -model.F_struc(:,2:end);
    prob.buc = full(model.F_struc(:,1));
    prob.blc = repmat(-inf,length(prob.buc),1);
else
    prob.a = sparse(ones(1,length(model.c))); % Dummy constraint
    prob.buc = inf;
    prob.blc = -inf;
end
if isempty(model.lb)
    prob.blx = repmat(-inf,1,length(model.c));
else
    prob.blx = model.lb;
end
if isempty(model.ub)
    prob.bux = repmat(inf,1,length(model.c));
else
    prob.bux = model.ub;
end

if model.K.f>0
    prob.blc(1:model.K.f) = prob.buc(1:model.K.f);
end

[prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(2*model.Q)));

if model.K.q(1)>0 || model.K.e > 0
    % Mosek crashes with empty lists, so only create if needed
    prob.cones.type = [];
    prob.cones.sub = [];
    prob.cones.subptr = [];
end

if model.K.q(1)>0
    % Append new variables to represent SOCP cones
    nof_new = sum(model.K.q);
    nof_original = size(prob.a,2);
    extendBasis = [spalloc(model.K.f,nof_new,0);
        spalloc(model.K.l,nof_new,0);
        speye(nof_new);
        spalloc(3*model.K.e,nof_new,0);
        spalloc(sum(model.K.s.^2),nof_new,0)];
    
    prob.a = [prob.a extendBasis];
    % Change the SOCP rows to be equalities for the slacks
    socpRows = 1+model.K.f+model.K.l:model.K.f+model.K.l+sum(model.K.q);
    prob.blc(socpRows) = prob.buc(socpRows);
    prob.c = [prob.c;zeros(nof_new,1)];
    % And now say that the new variables are in the SOCP cone
    top = nof_original;
    for i = 1:length(model.K.q)
        prob.cones.type = [prob.cones.type 0];
        prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+1];
        prob.cones.sub = [prob.cones.sub top+1:top+model.K.q(i)];
        prob.blx(top+1:top+model.K.q(i)) = -inf;
        prob.bux(top+1:top+model.K.q(i)) = inf;
        top = top + model.K.q(i);
    end
end

if nnz(model.K.e) > 0
    m = model.K.e;
    nof_new = m*3;
    nof_original = size(prob.a,2);
    extendedBasis = [spalloc(model.K.f + model.K.l + sum(model.K.q),nof_new,0);
        speye(m*3);
        spalloc(sum(model.K.s.^2),nof_new,0)];
    prob.a = [prob.a extendedBasis];
    prob.c = [prob.c;zeros(nof_new,1)];
    expRows = 1+model.K.f+model.K.l+sum(model.K.q):model.K.f+model.K.l+sum(model.K.q)+3*model.K.e;
    prob.blc(expRows) = prob.buc(expRows);
    prob.bux = [prob.bux;inf(nof_new,1)];
    prob.blx = [prob.blx;-inf(nof_new,1)];
    for i = 1:model.K.e
        prob.cones.type = [prob.cones.type 2];
        prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+1];
        prob.cones.sub    = [prob.cones.sub nof_original+3, nof_original+2, nof_original+1]; % YALMIPs exponential cone order different compare to mosek
        nof_original = nof_original + 3;
    end
end

if model.K.s(1)>0
    
    sdpRows = 1+model.K.f+model.K.l+3*model.K.e+sum(model.K.q):model.K.f+model.K.l+sum(model.K.q)+3*model.K.e+sum(model.K.s.^2);
    prob.blc(sdpRows) = prob.buc(sdpRows);
    
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subl = [];
    prob.bara.subk = [];
    prob.bardim  = model.K.s;
    prob.bara.val = [];
    
    prob.barc.subj = [];
    prob.barc.subk = [];
    prob.barc.subl = [];
    prob.barc.val = [];
    top = model.K.f + model.K.l + sum(model.K.q) + 3*model.K.e;
    A = prob.a(top+1:end,:);
    blc = prob.blc(top+1:end);
    buc = prob.buc(top+1:end);
    prob.a = prob.a(1:top,:);
    prob.buc = prob.buc(1:top,:);
    prob.blc = prob.blc(1:top,:);
    topA = 0;
    for i = 1:length(model.K.s)
        Z = (ones(model.K.s(i)));
        ii = find(Z);
        [kk,ll] = find(tril(Z));
        m = length(kk);
        jj = find(tril(Z));
        prob.a = [prob.a;A(topA + jj,:)];
        prob.buc = [prob.buc;buc(topA + jj)];
        prob.blc = [prob.blc;blc(topA + jj)];
        prob.bara.subi = [prob.bara.subi top + (1:m)];
        prob.bara.subj = [ prob.bara.subj  i*ones(1,m)];
        prob.bara.subk = [prob.bara.subk kk(:)'];
        prob.bara.subl = [ prob.bara.subl ll(:)'];
        temp = ones(1,m);
        temp(kk~=ll) = temp(kk~=ll)/2;
        prob.bara.val = [prob.bara.val temp];
        top = top + m;
        topA = topA + model.K.s(i)^2;
    end
end

if ~isempty(model.integer_variables)
    prob.ints.sub = model.integer_variables;
end

param = model.options.mosek;

% if ~isempty(model.x0) && model.K.s(1)==0 && model.K.q(1)==0
%     if model.options.usex0
% %         prob.sol.int.xx = zeros(max([length(model.Q) size(prob.a,2)]),1);
% %         prob.sol.int.xx(model.integer_variables) = model.x0(model.integer_variables);
% %         evalc('[r,res] = mosekopt (''symbcon'')');
% %         sc = res.symbcon ;
% %         param.MSK_IPAR_MIO_CONSTRUCT_SOL = sc.MSK_ON;
%     end
% end

% Debug?
if model.options.savedebug
    save mosekdebug prob param
end

if model.options.mosektaskfile
    mosekopt(sprintf('min write(%s) echo(0)', model.options.mosektaskfile), prob, param);
end

% Call MOSEK
showprogress('Calling MOSEK',model.options.showprogress);
if model.options.verbose == 0
    solvertime = tic;
    [r,res] = mosekopt('minimize echo(0)',prob,param);
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [r,res] = mosekopt('minimize',prob,param);
    solvertime = toc(solvertime);
end

if (r == 1010) || (r == 1011) | (r==1001)
    problem = -5;
    x = [];
    D_struc = [];
elseif r == 1200
    problem = 7;
    x = [];
    D_struc = [];
elseif r == 1295
    problem = -4;
    x = [];
    D_struc = [];
elseif r == 3100
    problem = 4;
    x = [];
    D_struc = [];
else
    % Recover solutions
    try
        sol = res.sol;
    catch
        error('Failure in Mosek recovery')
    end
    
    if isempty(model.integer_variables)
        x = sol.itr.xx(1:length(model.c)); % Might have added new ones
        D_struc = (sol.itr.suc-sol.itr.slc);
        error_message = sol.itr.prosta;
    else
        try
            x = sol.int.xx(1:length(model.c)); % Might have added new ones
            D_struc = [];
            error_message = sol.int.prosta;
        catch
            x = [];
            error_message = 'crash';
            D_struc = [];
        end
    end
    
    switch error_message
        case {'PRIMAL_AND_DUAL_FEASIBLE','PRIMAL_FEASIBLE'}
            problem = 0;
        case 'PRIMAL_INFEASIBLE'
            problem = 1;
        case 'DUAL_INFEASIBLE'
            problem = 2;
        case 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'
            problem = 12;
        otherwise
            problem = -1;
    end
end

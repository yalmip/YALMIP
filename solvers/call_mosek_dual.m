function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_dual(model)

% Convert if the caller is bnb or bmibnb which might have appended bounds
% Sure, we could setup model with bounds, but... 
[model.F_struc,model.K] = addStructureBounds(model.F_struc,model.K,model.ub,model.lb);

param = model.options.mosek;

prob.c = model.F_struc(1:model.K.f+model.K.l+sum(model.K.q)+3*model.K.e,1);
prob.a = -model.F_struc(1:model.K.f+model.K.l+sum(model.K.q)+3*model.K.e,2:end)';
prob.blc = -model.c;
prob.buc = -model.c;
prob.blx = -inf(size(prob.a,2),1);
prob.bux = inf(size(prob.a,2),1);
top = model.K.f+model.K.l;
prob.blx(1+model.K.f:model.K.f+model.K.l) = 0;

if model.K.q(1)>0 || model.K.e > 0
    prob.cones.type = [];
    prob.cones.subptr = [];
    prob.cones.sub = [];
end

if model.K.q(1)>0
    nq = length(model.K.q);
    prob.cones.type = zeros(nq, 1);
    prob.cones.subptr = zeros(nq, 1);
    prob.cones.sub = zeros(sum(model.K.q), 1);
    top0 = top;
    for i = 1:length(model.K.q)
        prob.cones.subptr(i) = top - top0 + 1;
        prob.cones.sub(top-top0+1:top-top0+model.K.q(i)) = top+1:top+model.K.q(i);
        top = top + model.K.q(i);
    end
end

if model.K.e>0
    for i = 1:model.K.e
        prob.cones.type = [prob.cones.type 3];
        prob.cones.subptr = [prob.cones.subptr length(prob.cones.sub)+1];
        prob.cones.sub = [prob.cones.sub top+3 top+2 top+1];        
        top = top + 3;
    end
end

if model.K.s(1)>0
   prob = appendMosekSDPdata(model.F_struc,model.K,prob);
end

if model.options.savedebug
    ops = model.options;
    save mosekdebug prob param
end

[r,res,solvertime] = doCall(prob,param,model.options);

try
    x = res.sol.itr.y;
catch   
    if isequal(model.options.mosek.MSK_IPAR_OPTIMIZER,'MSK_OPTIMIZER_FREE_SIMPLEX')
        x = res.sol.bas.y;
    else
        x = nan(length(model.c),1);    
    end
end

if model.options.saveduals & ~isempty(x)
    try       
        D_struc_SDP = zeros(sum(model.K.s.^2),1);
        top = 1;
        dtop = 1;
        for i = 1:length(model.K.s)          
            n = model.K.s(i);
            I = find(tril(ones(n)));
            v = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
            D_struc_SDP(dtop + I - 1) = v;
            in = ceil(I/n);
            jn = mod(I-1,n)+1;
            D_struc_SDP(dtop + (jn-1)*n+in - 1) = v;
            top = top + n*(n+1)/2;
            dtop = dtop  + n^2;
        end
        D_struc = [res.sol.itr.xx;D_struc_SDP];
    catch
        D_struc = [];
    end
else
    D_struc = [];
end

problem = MosekYALMIPError(res);

function [res,sol,solvertime] = doCall(prob,param,options)

showprogress('Calling Mosek',options.showprogress);
if options.verbose == 0
    solvertime = tic;
    [res,sol] = mosekopt('minimize echo(0)',prob,param);    
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [res,sol] = mosekopt('minimize info',prob,param);
    solvertime = toc(solvertime);
end

function problem = MosekYALMIPError(res)

if res.rcode == 2001
    problem = 1;
    return
elseif res.rcode == 1305
    problem = -4;
    return
elseif res.rcode == 10007
    problem = 16;
    return
elseif res.rcode == 1400
    problem = 20;
    return   
elseif res.rcode == 1001
    problem = -11;
    return;
elseif res.rcode == 1008
    problem = -12;
    return;
end

try
    solinfo = res.sol.itr;
catch
    solinfo = res.sol.bas;
end

switch solinfo.prosta
    case 'PRIMAL_AND_DUAL_FEASIBLE'        
        problem = 0;
    case 'DUAL_INFEASIBLE'
        problem = 1;
    case 'PRIMAL_INFEASIBLE'
        problem = 2;
    case 'MSK_RES_TRM_USER_CALLBACK'
        problem = 16;
    case 'MSK_RES_TRM_STALL'
        problem = 4;
    case 'UNKNOWN'
        try
            if isequal(res.rcodestr,'MSK_RES_TRM_STALL')
                problem = 4;
            elseif isequal(res.rcodestr,'MSK_RES_OK')
                problem = 0;
            else
                problem = 11;
            end
        catch
            problem = 9;
        end
    otherwise
        problem = -1;
end
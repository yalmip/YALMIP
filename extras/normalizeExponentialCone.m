function [model,output] = normalizeExponentialCone(model)

% N.B: YALMIP Definition % x2*exp(x1/x2)  <= x3

% Temporarily extract the data in a more standard name min c'*y, b + Ay >= 0
output.problem = 0;
if ~isempty(model.evalMap)

    % Temporarily extract the data in a more standard name min c'*y, b + Ay >= 0    
    data.A = model.F_struc(:,2:end);    
    data.b = full(model.F_struc(:,1));
    data.c = full(model.c);
    cones = model.K;
    % Clear away the SDP cone to easier add new rows. Will be put back in
    % place later
    data.A(end-sum(model.K.s.^2)+1:end,:)=[];
    data.b(end-sum(model.K.s.^2)+1:end) = [];
    sdpData = model.F_struc(end-sum(model.K.s.^2)+1:end,:);
    
    % First check that we have only exponentials
    convexFunctions = [];
    concaveFunctions = [];
    vectorentropy = [];
    vectorkullback = [];
    vectorlogsumexp = [];
    for i = 1:length(model.evalMap)
        switch model.evalMap{i}.fcn
            case {'exp','pexp'}
                convexFunctions = [convexFunctions model.evalMap{i}.computes];
            case {'log','plog','slog'}
                concaveFunctions = [concaveFunctions model.evalMap{i}.computes];
            case 'entropy'
                concaveFunctions = [concaveFunctions model.evalMap{i}.computes];
                if length(model.evalMap{i}.variableIndex) > 1
                    vectorentropy = [vectorentropy i];
                end
            case 'kullbackleibler'
                convexFunctions = [convexFunctions model.evalMap{i}.computes];
                if length(model.evalMap{i}.variableIndex) > 2
                    vectorkullback = [vectorkullback i];
                end
            case 'logsumexp'
                convexFunctions = [convexFunctions model.evalMap{i}.computes];
                vectorlogsumexp = [vectorlogsumexp i];
            otherwise
                % Standard interface, return solver not applicable
                % This should not be able to happen as we check this
                % earlier
                output = createoutput([],[],[],-4,model.solver.tag,[],[],0);
                return
        end
    end
    % Check that all exp/log enter in a convex fashion
    if model.K.f > 0
        if nnz(data.A(1:model.K.f,convexFunctions))>0 || nnz(data.A(1:model.K.f,concaveFunctions))>0
            output = createoutput([],[],[],19,model.solver.tag,[],[],0);
            return
        end
    end
    % Check sign in objective on exp/log terms
    if any(data.c(convexFunctions) < 0) || any(data.c(concaveFunctions) > 0)
        output = createoutput([],[],[],19,model.solver.tag,[],[],0);
        return
    end
    % Check sign in elementwise inequalities on exp/log terms
    if model.K.l > 0
        if nnz(data.A(1+model.K.f:model.K.f+model.K.l,convexFunctions)>0) || nnz(data.A(1+model.K.f:model.K.f+model.K.l,concaveFunctions)<0)
            output = createoutput([],[],[],19,model.solver.tag,[],[],0);
            return
        end
    end
    % Check so there are no exp/log terms in quadratic or SDP cone
    if sum(model.K.q) + sum(model.K.s) > 0
        if nnz(data.A(1+model.K.f+model.K.l:end,convexFunctions))>0 || nnz(data.A(1+model.K.f+model.K.l:end,concaveFunctions))>0
            output = createoutput([],[],[],19,model.solver.tag,[],[],0);
            return
        end
    end
    % Scalarize entropy operator
    if ~isempty(vectorentropy)
        for i = vectorentropy(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            n = length(data.c);
            m = length(involved);
            data.A = [sparse(ones(1,m+1),[computes n+(1:m)],[-1 ones(1,m)],1,n+m);
                      data.A spalloc(size(data.A,1),m,0)];
            data.b = [0;data.b];
            data.c = [data.c;zeros(m,1)];
            cones.f = cones.f + 1;
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = involved(1);
            % and then we add some
            for j = 2:length(involved)
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = involved(j);
                model.evalMap{end}.computes = n+j;
            end
        end
    end
    
    % Scalarize kullbackleibler operator
    if ~isempty(vectorkullback)
        for i = vectorkullback(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            n = length(data.c);
            m = length(involved)/2;
            data.A = [sparse(ones(1,m+1),[computes n+(1:m)],[-1 ones(1,m)],1,n+m);
                      data.A spalloc(size(data.A,1),m,0)];
            data.b = [0;data.b];
            data.c = [data.c;zeros(m,1)];
            cones.f = cones.f + 1;
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = involved([1 m+1]);
            % and then we add some
            for j = 2:length(involved)/2
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = involved([j j+m]);
                model.evalMap{end}.computes = n+j;
            end
        end
    end
    
    % Scalarize logsumexp operator
    % write log(expx1+expx2..)<= z as exp(x1-z)+exp(x2-z)+...<=1
    if ~isempty(vectorlogsumexp)
        for i = vectorlogsumexp(:)'
            involved = model.evalMap{i}.variableIndex;
            computes = model.evalMap{i}.computes;
            % Orignal #variables
            n = length(data.c);
            % Number of terms in logsumexp
            m = length(involved);
            
            % exp(x1-z)+exp(x2-z)+...<=1
            data.A = [data.A(1:cones.f,:) spalloc(cones.f,m,0);
                      spalloc(1,n,0) -ones(1,m);
                      data.A(cones.f+1:end,:) spalloc(cones.l+sum(cones.q)+sum(cones.s),m,0)];
            data.b = [data.b(1:cones.f);
                      1;
                      data.b(cones.f+1:end)];
            data.c = [data.c;zeros(m,1)];
            cones.l = cones.l + 1;
            
            % The original strucuture now just relates to the first new term
            model.evalMap{i}.computes = n+1;
            model.evalMap{i}.variableIndex = [involved(1) computes];
            model.evalMap{i}.fcn = 'expdiff';
            % and then we add some new
            for j = 2:length(involved)
                model.evalMap{end+1} = model.evalMap{i};
                model.evalMap{end}.variableIndex = [involved(j) computes];
                model.evalMap{end}.computes = n+j;
            end
        end
    end
    
    % Describe all exponential cones
    m = length(model.evalMap);
    cones.e = model.K.e + m;       
    for i = 1:m
        switch model.evalMap{i}.fcn
            case 'exp'
                % 1*exp(xv/1) <= xc
                % y*exp(x/y)  <= z.  y new variable, xv the original variable, and xc the "computed"
                x = model.evalMap{i}.variableIndex;
                z = model.evalMap{i}.computes;
                data.A = [data.A;
                          -sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];
            case 'pexp'
                % xv(1)*exp(xv(2)/xv(1)) <= xc
                x = model.evalMap{i}.variableIndex(2);
                y = model.evalMap{i}.variableIndex(1);
                z = model.evalMap{i}.computes;
                data.A = [data.A;
                          -sparse([1;2;3],[x y z],[-1 -1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]];
            case 'log'
                % log(xv) >= xc i.e. xv >= exp(xc/1)*1
                z = model.evalMap{i}.variableIndex;
                x = model.evalMap{i}.computes;
                data.A = [data.A;
                          -sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];
            case 'plog'
                % xv(1)log(xv(2)/xv(1))>=xc i.e.
                % -xc >= -xv(1)log(xv(2)/xv(1))
                z = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex(1);
                x = model.evalMap{i}.variableIndex(2);
                data.A = [data.A;
                          -sparse([1;2;3],[x y z],[1 -1 1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]];
            case 'slog'
                % log(1+xv) >= xc i.e. (1+xv) >= exp(xc/1)*1
                z = model.evalMap{i}.variableIndex;
                x = model.evalMap{i}.computes;
                data.A = [data.A;
                          -sparse([1;3],[x z],-1,3,size(data.A,2))];
                data.b = [data.b;0;1;1];
            case 'entropy'
                % -xv*log(xv)>=xc i.e. 1 >= exp(xc/xv)*xv
                x = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex;
                data.A = [data.A;
                          -sparse([1;2],[x y],[-1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;1]];
            case 'kullbackleibler'
                % -xv1*log(xv1/xv2)>=-xc i.e. xv2 >= exp(-xc/xv1)*xv1
                x = model.evalMap{i}.computes;
                y = model.evalMap{i}.variableIndex(1);
                z = model.evalMap{i}.variableIndex(2);
                data.A = [data.A;
                          -sparse([1;2;3],[x y z],[1 -1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;0;0]];
            case 'expdiff'
                % exp(xv(1)-xv(2)) <= xc
                x = model.evalMap{i}.variableIndex;
                z = model.evalMap{i}.computes;
                data.A = [data.A;
                          -sparse([1;1;3],[x(1) x(2) z],[-1 1 -1],3,size(data.A,2))];
                data.b = [data.b;[0;1;0]];
            otherwise
        end
    end
    model.F_struc = [data.b data.A];
    model.c = data.c;
    model.K = cones;   
    % Put back the SDP cone in place
    if model.K.s(1)>0
        if size(sdpData,2) < size(model.F_struc,2)
            sdpData(end,size(model.F_struc,2)) = 0;
        end
        model.F_struc = [model.F_struc;sdpData];
    end
end
n = size(model.F_struc,2)-1;
if n > length(model.lb)
    missing = n - length(model.lb);
    model.lb = [model.lb;-inf(missing,1)];
end
if size(model.F_struc,2)-1 > length(model.ub)
    missing = n - length(model.ub);
    model.ub = [model.ub;inf(missing,1)];
end

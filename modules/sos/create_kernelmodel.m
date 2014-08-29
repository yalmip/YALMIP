function [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,sparsityPattern);

% To get the primal kernel representation, we simply use
% the built-in dualization module!
% First, write the problem in primal kernel format
traceobj = 0;
dotraceobj = options.sos.traceobj;
F = F_parametric;
for i = 1:length(Blockedb)


    sizematrixSOS = sqrt(size(Blockedb{i},2));
    for k = 1:sizematrixSOS
        for r = k:sizematrixSOS
            res{(k-1)*sizematrixSOS+r} = 0;
        end
    end

    for j = 1:length(BlockedA{i})
        n = sqrt(size(BlockedA{i}{j},2));
        BlockedQ{i}{j} = sdpvar(n*sizematrixSOS,n*sizematrixSOS);
        F = F + lmi(BlockedQ{i}{j});
        if sizematrixSOS>0
            % Matrix valued sum of sqaures
            % Loop over all elements
            starttop = 1;
            for k = 1:sizematrixSOS
                startleft = 1;
                for r = 1:sizematrixSOS
                    if k<=r
                        Qkr = BlockedQ{i}{j}(starttop:starttop+n-1,startleft:startleft+n-1);
                        res{(k-1)*sizematrixSOS+r} = res{(k-1)*sizematrixSOS+r} + BlockedA{i}{j}*reshape(Qkr,n^2,1);
                    end
                    startleft = startleft + n;
                end
                starttop = starttop + n;
            end
        else
            % Standard case
            res{1} = res{1} + BlockedA{i}{j}*reshape(BlockedQ{i}{j},n^2,1);
        end
        if dotraceobj
            traceobj = traceobj + trace(BlockedQ{i}{j});
        end
    end
    for k = 1:sizematrixSOS
        for r = k:sizematrixSOS
            F = F + (res{(k-1)*sizematrixSOS+r} == Blockedb{i}(:,(k-1)*sizematrixSOS+r));
        end
    end
end

% % Constrain elements according to a desired sparsity
if ~isempty(sparsityPattern)
    res = [];
    for i = 1:length(BlockedQ)
        for j = 1:length(BlockedQ{i})
            if ~isempty(sparsityPattern{i}{j})
                H = spalloc(length(BlockedQ{i}{j}),length(BlockedQ{i}{j}),length(sparsityPattern{i}{j}));
                H(sparsityPattern{i}{j}) = 1;
                k = find(triu(H));
                res = [res;BlockedQ{i}{j}(k)];
            end
        end
    end
    F = F + (0 == res);
end

% And get the primal model of this
if isempty(parobj)
    if options.sos.traceobj
        [F,obj,Primal_matrices,Free_variables] = dualize(F,traceobj,1,options.sos.extlp);
    else
        [F,obj,Primal_matrices,Free_variables] = dualize(F,[],1,options.sos.extlp);
    end
else
    [F,obj,Primal_matrices,Free_variables] = dualize(F,parobj,1,options.sos.extlp);
end
% In dual mode, we maximize
obj = -obj;

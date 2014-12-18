function [F,obj,BlockedQ,sol] = create_imagemodel(BlockedA,Blockedb,F_parametric,parobj,options);


% *********************************
% IMAGE REPRESENTATION
% Needed for nonlinearly parameterized problems
% More efficient in some cases
% *********************************
g = [];
vecQ = [];
sol.problem = 0;
for i = 1:length(BlockedA)
    q = [];
    A = [];
    for j = 1:length(BlockedA{i})
        n = sqrt(size(BlockedA{i}{j},2));
        BlockedQ{i}{j} = sdpvar(n,n);
        q = [q;reshape(BlockedQ{i}{j},n^2,1)];
        A = [A BlockedA{i}{j}];
    end
    vecQ = [vecQ;q];
    g = [g;Blockedb{i}-A*q];
end

g_vars = getvariables(g);
q_vars = getvariables(vecQ);
x_vars = setdiff(g_vars,q_vars);

base = getbase(g);
if isempty(x_vars)
    A = base(:,1);base = base(:,2:end);
    B = (base(:,ismember(g_vars,q_vars)));
    Bnull = sparsenull(B);
    t = sdpvar(size(Bnull,2),1);
    imQ = -B\A+Bnull*t;
else
    A = base(:,1);base = base(:,2:end);
    C = base(:,ismember(g_vars,x_vars));
    B = (base(:,ismember(g_vars,q_vars)));
    [Bnull,Q1,R1] = sparsenull(B);
    [ii,jj,kk] = find(Bnull);
    %Bnull(abs(Bnull) < 1e-12) = 0;
    small = find(abs(kk)<1e-12);
    ii(small)=[];
    jj(small)=[];
    kk(small)=[];
    Bnull = sparse(ii,jj,kk,size(Bnull,1),size(Bnull,2));    
    t = sdpvar(size(Bnull,2),1);
    if options.sos.model == 2
        imQ = -Q1*(R1'\(A+C*recover(x_vars)))+Bnull*t;
    else
        H = deriveBasis(B);
        imQ = -B\(A+C*recover(x_vars))+H*t;        
    end
end
notUsed = find(sum(abs(B),2)==0);
if ~isempty(notUsed)
    ff=g(notUsed);
    if isa(ff,'double')
        if nnz(ff)>0
            sol.yalmiptime = 0; % FIX
            sol.solvertime = 0;
            sol.problem = 2;
            sol.info = yalmiperror(1,'YALMIP');
            F = [];
            obj = [];
            if options.verbose > 0
                disp(' ');
                disp('-> During construction of data, I encountered a situation')
                disp('-> situation that tells me that the problem is trivially infeasible!')
                disp('-> Have you forgotten to define some parametric variables?,')
                disp('-> or perhaps you have a parametric problem where the highest')
                disp('-> power in some of the independent variables is odd for any choice')
                disp('-> of parametric variables, such as x^8+x^7+t*y^3')
                disp('-> Anyway, take a good look at your model again...');                
            end
            return
            %            error('You seem to have a strange model. Have you forgotten to define some parametric variable?');
        end
    else
       F_parametric = F_parametric + (g(notUsed)==0);
    end
end
F_sos = ([]);
obj = 0;
for i = 1:length(BlockedQ)
    for j = 1:size(BlockedQ{i},2)
        Q_old = BlockedQ{i}{j};
        Q_old_vars = getvariables(Q_old);
        Q_old_base = getbase(Q_old);
        in_this = [];
        for k = 1:length(Q_old_vars)
            in_this = [in_this;find(Q_old_vars(k)==q_vars)];
        end
        Q_new = Q_old_base(:,1) + Q_old_base(:,2:end)*imQ(in_this);
        Q_new = reshape(Q_new,length(BlockedQ{i}{j}),length(BlockedQ{i}{j}));
        obj = obj+trace(Q_new);
        if ~isa(Q_new,'double')
            switch options.sos.model
                case {2,4}
                    F_sos = F_sos + [Q_new>=0];
                case 5                    
                    F_sos = F_sos + [dd(Q_new)];
                case 6
                    F_sos = F_sos + [sdd(Q_new)];                                             
                otherwise
            end
        elseif min(eig(Q_new))<-1e-8
            sol.yalmiptime = 0; % FIX
            sol.solvertime = 0;
            sol.problem = 2;
            sol.info = yalmiperror(1,'YALMIP');
            F = [];
            obj = [];
            error('Problem is trivially infeasible. After block-diagonalizing, I found constant negative definite blocks!');
            return
        end
        BlockedQ{i}{j} = Q_new;
    end
end

F = F_parametric + F_sos;

if isa(obj,'double') | (options.sos.traceobj == 0)
    obj = [];
end

if ~isempty(parobj)
    obj = parobj;
end
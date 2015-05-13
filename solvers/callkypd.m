function diagnostic = callkypd(F,h,options)

%F = constraint2kyp(F);
kyps = is(F,'kyp');
if all(kyps)
    kypConstraints = F;
    otherConstraints = lmi;
else
    kypConstraints   = F(find(kyps));
    otherConstraints = F(find(~kyps));
end

% Trivial case
if length(kypConstraints) == 0
    if ~isempty(options)
        options.solver = '';
    else
        options = sdpsettings;
    end
    diagnostic = solvesdp(F,h,options)
    return;
end

p_variables = [];
x_variables = [];
for i = 1:length(kypConstraints)
    [A{i},B{i},P{i},M{i},negated{i}] = extractkyp(sdpvar(kypConstraints(i)));
    if ~isa(P{i},'sdpvar') | ~isempty(intersect(p_variables,getvariables(P{i})))
        diagnostic.solvertime = 0;
        diagnostic.problem = -4;
        diagnostic.info = yalmiperror(-4);   
        return;
    end
    P_base = getbase(P{i});
    [n,m] = size(P{i});
    if (n~=m) | (nnz(P_base)~=n*n) | (sum(sum(P_base))~=n*n)
        diagnostic.solvertime = 0;
        diagnostic.problem = -4;
        diagnostic.info = yalmiperror(-4);   
        return;        
    end
    p_variables = [p_variables getvariables(P{i})];
    x_variables = [x_variables getvariables(M{i})];
end

if intersect(p_variables,getvariables(h))
    diagnostic.solvertime = 0;
    diagnostic.problem = -4;
    diagnostic.info = yalmiperror(-4);   
    return;
end

start = length(A);k = 1;
for i = 1:length(otherConstraints)
    if is(otherConstraints(i),'lmi')
        A{k+start}=[];
        B{k+start}=[];
        P{k+start}=[];
        M{k+start} = sdpvar(otherConstraints(i));
        x_variables = [x_variables getvariables(M{k+start})];
        k = k + 1;
    elseif is(otherConstraints(i),'elementwise')
        X = sdpvar(otherConstraints(i));
        for ii = 1:size(X,1)
            for jj = 1:size(X,2)
                A{k+start}=[];
                B{k+start}=[];
                P{k+start}=[];
                M{k+start} = X(ii,jj);
                x_variables = [x_variables getvariables(X(ii,jj))];
                k = k + 1;
            end
        end
    end
end

% Objective function variables
x_variables = unique([x_variables getvariables(h)]);

% Are P and x independent
if ~isempty(intersect(p_variables,x_variables))
    error
end

% Generate the M-matrices
for i = 1:length(M)
    m_variables = getvariables(M{i});
    n = length(M{i});
    M0{i} = getbasematrix(M{i},0);
    for j = 1:length(x_variables)
        Mi{i,j} = getbasematrix(M{i},x_variables(j));
    end
end

c = zeros(length(x_variables),1);
if ~isempty(h)
    for i = 1:length(x_variables)
        c(i) = getbasematrix(h,x_variables(i));
    end
end

for i = 1:length(M)
    matrixinfo.n(i) = length(A{i});
    matrixinfo.nm(i) = length(M{i});
end

matrixinfo.K  = length(x_variables);
matrixinfo.c  = c;
matrixinfo.N  = length(A);
matrixinfo.A  = A;
matrixinfo.B  = B;
matrixinfo.M0 = M0;
matrixinfo.M  = Mi;
matrixinfo.A  = A;

try
    solvertime = tic; 
    [u,Popt,xopt,Zopt] = kypd_solver(matrixinfo,options);
    solvertime = toc(solvertime);
    
    % SAVE PRIMALS
    setsdpvar(recover(x_variables),xopt);
    for i = 1:length(kypConstraints)
        if negated{i}
            setsdpvar(P{i},-Popt{i});
        else
            setsdpvar(P{i},Popt{i});
        end
    end
    
    % SAVE DUALS
    lmi_index = [];
    j = 1;
    for i = 1:length(kypConstraints)
        %W = struct(kypConstraints(i));lmiid = W.LMIid;
        lmiid = getlmiid(kypConstraints(i));        
        lmi_index = [lmi_index;lmiid];
        duals{j}=Zopt{j};j = j+1;
    end
    for i = 1:length(otherConstraints)
        %W = struct(otherConstraints(i));lmiid = W.LMIid;
        lmiid = getlmiid(otherConstraints(i));        
        lmi_index = [lmi_index;lmiid];
        duals{j}=Zopt{j};j = j+1;
    end       
    yalmip('setdual',lmi_index,duals);
    
    diagnostic.solvertime = solvertime;
    diagnostic.info = yalmiperror(0,'KYPD');
    diagnostic.problem = 0;
catch
    solvertime = etime(clock,solvertime);    
    diagnostic.solvertime = solvertime;
    diagnostic.info = yalmiperror(9,'KYPD');
    diagnostic.problem = 9;
end


function Fout = constraint2kyp(F)
Fout = [];
for i = 1:length(F)
    if is(F(i),'lmi')
        [l,m,r] = factors(sdpvar(F(i)));
        [constants,general,singles,pairs] = classifyfactors(l,m,r);
        [constants,general,singles,pairs] = compressfactors2(constants,general,singles,pairs);    
    else
        Fout = Foput + F(i);
    end
end











function [F_struc,K,KCut,schur_funs,schur_data,schur_variables] = lmi2sedumistruct(F)
%lmi2sedumistruct   Internal function: Converts LMI to format needed in SeDuMi

nvars = yalmip('nvars'); %Needed lot'sa times...

% We first browse to see what we have got and the
% dimension of F_struc (massive speed improvement)
F = flatten(F);
type_of_constraint = zeros(length(F.LMIid),1);%zeros(size(F.clauses,2),1);
any_cuts = 0;
for i = 1:length(F.LMIid)%size(F.clauses,2)
    type_of_constraint(i) = F.clauses{i}.type;
    if F.clauses{i}.cut
        any_cuts = 1;
    end
end

F_struc = [];

schur_data = [];
schur_funs = [];
schur_variables = [];

sdp_con = find(type_of_constraint == 1 | type_of_constraint == 9 | type_of_constraint == 40);
sdpstack_con = find(type_of_constraint == 57);
lin_con = find(type_of_constraint == 2 | type_of_constraint == 12);
equ_con = find(type_of_constraint == 3);
qdr_con = find(type_of_constraint == 4);
mqdr_con = find(type_of_constraint == 54);
rlo_con = find(type_of_constraint == 5);
pow_con = find(type_of_constraint == 20);
sos2_con = find(type_of_constraint == 50);
sos1_con = find(type_of_constraint == 51);
cmp_con = find(type_of_constraint == 55);

% SeDuMi struct
K.f = 0; % Linear equality
K.l = 0; % Linear inequality
K.c = 0; % Complementarity constraints
K.q = 0; % SOCP
K.r = 0; % Rotated SOCP (obsolete)
K.p = 0; % Power cone
K.s = 0; % SDP cone
K.rank = [];
K.dualrank = [];
K.scomplex = [];
K.xcomplex = [];

KCut.f = [];
KCut.l = [];
KCut.c = [];
KCut.q = [];
KCut.r = [];
KCut.p = [];
KCut.s = [];

top = 1;
localtop = 1;

% In the first part of the code, we will work with a transposed version of
% the vectorized constraints, i.e. constraints are added by adding columns,
% wheras the final output will be transposed

% Linear equality constraints
alljx = [];
allix = [];
allsx = [];
block = 0;
for i = 1:length(equ_con)
    constraints = equ_con(i);
    data = getbase(F.clauses{constraints}.data);
  %  [n,m] = size(F.clauses{constraints}.data);
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    if isreal(data)
        ntimesm = size(data,1);
        %ntimesm = n*m; %Just as well pre-calc
    else
        % Complex constraint, Expand to real and Imag
        ntimesm = 2*size(data,1);
        %ntimesm = 2*n*m; %Just as well pre-calc
        data = [real(data);imag(data)];
    end
    mapX = [1 1+lmi_variables];
    [ix,jx,sx] = find(data);    
   
    %F_structemp = sparse(mapX(jx),ix,sx,1+nvars,ntimesm);    
    %F_struc = [F_struc F_structemp];
    
    alljx = [alljx mapX(jx)];
    allix = [allix ix(:)'+block];block = block + ntimesm;
    allsx = [allsx sx(:)'];
    
    if F.clauses{constraints}.cut
        KCut.f = [KCut.f localtop:localtop+ntimesm-1];
    end
    
    localtop = localtop+ntimesm;
    top = top+ntimesm;
    K.f = K.f+ntimesm;
end
F_struc = sparse(alljx,allix,allsx,1+nvars,block);

% Linear inequality constraints
localtop = 1;

% Cuts are not dealt with correctly in the recurisve code. Use slower
% version. Test with bmibnb_qcqp5
if any_cuts
    for i = 1:length(lin_con)
        constraints = lin_con(i);
        data = getbase(F.clauses{constraints}.data);
        [n,m] = size(F.clauses{constraints}.data);
        
        % Which variables are needed in this constraint
        lmi_variables = getvariables(F.clauses{constraints}.data);
        
        % Convert to real problem
        if isreal(data)
            ntimesm = n*m; %Just as well pre-calc
        else
            % Complex constraint, Expand to real and Imag
            ntimesm = 2*n*m; %Just as well pre-calc
            data = [real(data);imag(data)];
        end
        
        % Add numerical data to complete problem setup
        mapX = [1 1+lmi_variables];
        [ix,jx,sx] = find(data);
        F_structemp = sparse(mapX(jx),ix,sx,1+nvars,ntimesm);
        F_struc = [F_struc F_structemp];
        
        if F.clauses{constraints}.cut
            KCut.l = [KCut.l localtop:localtop+ntimesm-1];
        end
        
        localtop = localtop+ntimesm;
        top = top+ntimesm;
        K.l = K.l+ntimesm;
    end
else
    if length(lin_con)>0
        [F_struc,K,KCut] = recursive_lp_fix(F,F_struc,K,KCut,lin_con,nvars,4,1);
    end
end

for i = 1:length(cmp_con)
    constraints = cmp_con(i);
    
    [n,m] = size(F.clauses{constraints}.data);
    ntimesm = n*m; %Just as well pre-calc
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    
    % We allocate the structure blockwise...
    F_structemp  = spalloc(1+nvars,ntimesm,0);
    % Add these rows only
    F_structemp([1 1+lmi_variables(:)'],:)= getbase(F.clauses{constraints}.data).';
    
    % ...and add them together (efficient for large structures)
    F_struc = [F_struc F_structemp];
    
    top = top+ntimesm;
    K.c(i) = n;
end

qdr_con = union(qdr_con,mqdr_con);
if length(qdr_con) > 0
    [F_struc,K,KCut] = recursive_socp_fix(F,F_struc,K,KCut,qdr_con,nvars,8,1);
end
% % 
% if length(mqdr_con) > 0
%    [F_struc,K,KCut] = recursive_msocp_fix(F,F_struc,K,KCut,mqdr_con,nvars,inf,1);
% end


% Rotated Lorentz cone constraints
for i = 1:length(rlo_con)
    constraints = rlo_con(i);
    
    [n,m] = size(F.clauses{constraints}.data);
    ntimesm = n*m; %Just as well pre-calc
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    
    % We allocate the structure blockwise...
    F_structemp  = spalloc(1+nvars,ntimesm,0);
    % Add these rows only
    F_structemp([1 1+lmi_variables(:)'],:)= getbase(F.clauses{constraints}.data).';
    
    % ...and add them together (efficient for large structures)
    F_struc = [F_struc F_structemp];
    
    top = top+ntimesm;
    K.r(i) = n;
end

% Power cone constraints
for i = 1:length(pow_con)
    constraints = pow_con(i);
    
    [n,m] = size(F.clauses{constraints}.data);
    ntimesm = n*m; %Just as well pre-calc
    
    % Should always have size 4
    if n~=4
        error('Power cone constraint has strange dimension')
    end
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    
    % We allocate the structure blockwise...
    F_structemp  = spalloc(1+nvars,ntimesm,0);
    % Add these rows only
    F_structemp([1 1+lmi_variables(:)'],:)= getbase(F.clauses{constraints}.data).';
    
    alpha = F_structemp(1,end);
    F_structemp(:,end)=[];
    % ...and add them together (efficient for large structures)
    F_struc = [F_struc F_structemp];
    
    top = top+ntimesm;
    K.p(i) = alpha;
end



% Semidefinite  constraints
% We append the recursively in order to speed up construction
% of problems with a lot of medium size SDPs
any_schur = 0;
for j = sdp_con(:)'
    if ~isempty(F.clauses{j}.schurfun)
        any_schur = 1;
        break
    end
end
if any_schur
    [F_struc,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,F_struc,K,KCut,schur_funs,schur_data,schur_variables,sdp_con,nvars,inf,1);
else
    [F_struc,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,F_struc,K,KCut,schur_funs,schur_data,schur_variables,sdp_con,nvars,8,1);
end
% Now go back to YALMIP default orientation constraints x variables
F_struc = F_struc.';

% Now fix things for the rank constraint
% This is currently a hack...
% Should not be in this file
[rank_variables,dual_rank_variables] = yalmip('rankvariables');
if ~isempty(rank_variables)
    used_in = find(sum(abs(F_struc(:,1+rank_variables)),2));
    if ~isempty(used_in)
        if used_in >=1+K.f & used_in < 1+K.l+K.f
            for i = 1:length(used_in)
                [ii,jj,kk] = find(F_struc(used_in(i),:));
                if length(ii)==2 & kk(2)<1
                    r = floor(kk(1));
                    var = jj(2)-1;
                    extstruct = yalmip('extstruct',var);
                    X = extstruct.arg{1};
                    if issymmetric(X)
                        F_structemp = sedumize(X,nvars);
                    else
                        error('Only symmetric matrices can be rank constrained.')
                    end
                    F_struc = [F_struc;F_structemp];
                    if isequal(K.s,0)
                        K.s(1,1) = size(extstruct.arg{1},1);
                    else
                        K.s(1,end+1) = size(extstruct.arg{1},1);
                    end
                    K.rank(1,end+1) = min(r,K.s(end));
                else
                    error('This rank constraint is not supported (only supports rank(X) < r)')
                end
            end
            % Remove the nonlinear operator constraints
            
            F_struc(used_in,:) = [];
            K.l = K.l - length(used_in);
        else
            error('You have added a rank constraint on an equality constraint, or a scalar expression?!')
        end
    end
end
if ~isempty(sos2_con) | ~isempty(sos1_con)
    K.sos.type  = [];
    K.sos.variables  = {};
    K.sos.weight  = {};
    for i = sos2_con(:)'
        K.sos.type = [K.sos.type '2'];
        K.sos.variables{end+1} = getvariables(F.clauses{i}.data);
        K.sos.variables{end} = K.sos.variables{end}(:);
        temp = struct(F.clauses{i}.data);
        K.sos.weight{end+1} = temp.extra.sosweights;
    end
    for i = sos1_con(:)'
        K.sos.type = [K.sos.type '1'];
        K.sos.variables{end+1} = getvariables(F.clauses{i}.data);
        K.sos.variables{end} = K.sos.variables{end}(:);
        temp = struct(F.clauses{i}.data);
        K.sos.weight{end+1} = temp.extra.sosweights;
    end
else
    K.sos.type = [];
    K.sos.variables = [];
    K.sos.weight = [];
end
if ~isempty(dual_rank_variables)
    used_in = find(sum(abs(F_struc(:,1+dual_rank_variables)),2));
    if ~isempty(used_in)
        if used_in >=1+K.f & used_in < 1+K.l+K.f
            for i = 1:length(used_in)
                [ii,jj,kk] = find(F_struc(used_in(i),:));
                if length(ii)==2 & kk(2)<1
                    r = floor(kk(1));
                    var = jj(2)-1;
                    extstruct = yalmip('extstruct',var);
                    X = extstruct.arg{1};
                    id = getlmiid(X);
                    inlist=getlmiid(F);
                    index=find(id==inlist);
                    if ~isempty(index)
                        K.rank(1,index) = min(r,K.s(index));
                    end
                else
                    error('This rank constraint is not supported (only supports rank(X) < r)')
                end
            end
            % Remove the nonlinear operator constraints
            F_struc(used_in,:) = [];
            K.l = K.l - length(used_in);
        else
            error('You have added a rank constraint on an equality constraint, or a scalar expression?!')
        end
    end
end

function F_structemp = sedumize(Fi,nvars)
Fibase = getbase(Fi);
[n,m] = size(Fi);
ntimesm = n*m;
lmi_variables = getvariables(Fi);
[ix,jx,sx] = find(Fibase);
mapX = [1 1+lmi_variables];
F_structemp = sparse(ix,mapX(jx),sx,ntimesm,1+nvars);

function [F_struc,K,KCut] = recursive_lp_fix(F,F_struc,K,KCut,lp_con,nvars,maxnlp,startindex)

% Check if we should recurse
if length(lp_con)>=2*maxnlp
    % recursing costs, so do 4 in one step
    ind = 1+ceil(length(lp_con)*(0:0.25:1));
    [F_struc1,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(ind(1):ind(2)-1),nvars,maxnlp,startindex+ind(1)-1);
    [F_struc2,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(ind(2):ind(3)-1),nvars,maxnlp,startindex+ind(2)-1);
    [F_struc3,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(ind(3):ind(4)-1),nvars,maxnlp,startindex+ind(3)-1);
    [F_struc4,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(ind(4):ind(5)-1),nvars,maxnlp,startindex+ind(4)-1);
    F_struc = [F_struc F_struc1 F_struc2 F_struc3 F_struc4];
    return
elseif length(lp_con)>=maxnlp
    mid = ceil(length(lp_con)/2);
    [F_struc1,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(1:mid),nvars,maxnlp,startindex);
    [F_struc2,K,KCut] = recursive_lp_fix(F,[],K,KCut,lp_con(mid+1:end),nvars,maxnlp,startindex+mid);
    F_struc = [F_struc F_struc1 F_struc2];
    return
end

oldF_struc = F_struc;
F_struc = [];
for i = 1:length(lp_con)
    constraints = lp_con(i);
    Fi = F.clauses{constraints}.data;
    Fibase = getbase(Fi);
  %  [n,m] = size(Fi);
    
    % Convert to real problem
    if isreal(Fibase)
        ntimesm = size(Fibase,1);
        %ntimesm = n*m; %Just as well pre-calc
    else
        % Complex constraint, Expand to real and Imag
        ntimesm = 2*size(Fibase,1);
        %ntimesm = 2*n*m; %Just as well pre-calc
        Fibase = [real(Fibase);imag(Fibase)];
    end
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(Fi);
    mapX = [1 1+lmi_variables];
    
  %  simpleMap = all(mapX==1:length(mapX));
    simpleMap =  0;%all(diff(lmi_variables)==1);
    % highly optimized concatenation...
    if size(Fibase) == [ntimesm 1+nvars] & simpleMap     
        F_struc = [F_struc Fibase'];
    elseif simpleMap
        vStart = lmi_variables(1);
        vEnd = lmi_variables(end);
        if vStart == 1
            F_struc = [F_struc [Fibase spalloc(n,nvars-vEnd,0)]'];
        else
            F_struc = [F_struc [Fibase(:,1) spalloc(n,vStart-1,0) Fibase(:,2:end) spalloc(n,nvars-vEnd,0)]'];
        end
    else
        [ix,jx,sx] = find(Fibase);
        F_structemp = sparse(mapX(jx),ix,sx,1+nvars,ntimesm);
        F_struc = [F_struc F_structemp];
    end
    
    if F.clauses{constraints}.cut
        KCut.l = [KCut.l i+startindex-1:i+startindex-1+n];
    end
    
    K.l(i+startindex-1) = ntimesm;
end
K.l = sum(K.l);
F_struc = [oldF_struc F_struc];


function [F_struc,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,F_struc,K,KCut,schur_funs,schur_data,schur_variables,sdp_con,nvars,maxnsdp,startindex)

if isempty(sdp_con)
    return
    % Check if we should recurse
elseif length(sdp_con)>2*maxnsdp
    % recursing costs, so do 4 in one step
    ind = 1+ceil(length(sdp_con)*(0:0.25:1));
    [F_struc1,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(ind(1):ind(2)-1),nvars,maxnsdp,startindex+ind(1)-1);
    [F_struc2,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(ind(2):ind(3)-1),nvars,maxnsdp,startindex+ind(2)-1);
    [F_struc3,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(ind(3):ind(4)-1),nvars,maxnsdp,startindex+ind(3)-1);
    [F_struc4,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(ind(4):ind(5)-1),nvars,maxnsdp,startindex+ind(4)-1);
    F_struc = [F_struc F_struc1 F_struc2 F_struc3 F_struc4];
    return
elseif length(sdp_con)>maxnsdp
    mid = ceil(length(sdp_con)/2);
    [F_struc1,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(1:mid),nvars,maxnsdp,startindex);
    [F_struc2,K,KCut,schur_funs,schur_data,schur_variables] = recursive_sdp_fix(F,[],K,KCut,schur_funs,schur_data,schur_variables,sdp_con(mid+1:end),nvars,maxnsdp,startindex+mid);
    F_struc = [F_struc F_struc1 F_struc2];
    return
end

oldF_struc = F_struc;
F_struc = [];
for i = 1:length(sdp_con)
    constraints = sdp_con(i);
    
    % Simple data
    Fi = F.clauses{constraints};
    X = Fi.data;
    lmi_variables = getvariables(X);
    [n,m] = size(X);
    ntimesm = n*m; %Just as well pre-calc
    
    if is(X,'gkyp')
        ss = struct(X);
        
        nn = size(F.clauses{1}.data,1);
        bb = getbase(ss.extra.M);
        Mbase = [];
        for ib = 1:size(bb,2)-1
            Mbase{ib} = (reshape(bb(:,ib+1),nn,nn));
        end
        for ip = 1:length(ss.extra.P)
            Pvars{ip} = getvariables(ss.extra.P{ip});
        end
        schur_data{i,1} = {ss.extra.K, ss.extra.Phi,Mbase, ss.extra.negated, getvariables(ss.extra.M),Pvars};
        schur_funs{i,1} = 'HKM_schur_GKYP';
        schur_variables{i,1} = lmi_variables;
        
    elseif ~isempty(Fi.schurfun)
        schur_data{i,1} = Fi.schurdata;
        schur_funs{i,1} = Fi.schurfun;
        schur_variables{i,1} = lmi_variables;
    end
    
    % get numerics
    Fibase = getbase(X);
    % now delete old data to save memory
    % F.clauses{constraints}.data=[];
    
    % Which variables are needed in this constraint
    %lmi_variables = getvariables(Fi);
    if length(lmi_variables) == nvars
        % No remap needed
        F_structemp =  Fibase';
    else
        mapX = [1 1+lmi_variables];
        [ix,jx,sx] = find(Fibase);
        clear Fibase;
        % Seems to be faster to transpose generation
        F_structemp = sparse(ix,mapX(jx),sx,ntimesm,1+nvars)';
        clear jx ix sx
    end
    F_struc = [F_struc F_structemp];
    
    if Fi.cut
        KCut.s = [KCut.s i+startindex-1];
    end
    K.s(i+startindex-1) = n;
    K.rank(i+startindex-1) = n;
    K.dualrank(i+startindex-1) = n;
    % Check for a complex structure
    if ~isreal(F_structemp)
        K.scomplex = [K.scomplex i+startindex-1];
    end
    clear F_structemp
    
end
F_struc = [oldF_struc F_struc];





function [F_struc,K,KCut] = recursive_socp_fix(F,F_struc,K,KCut,qdr_con,nvars,maxnsocp,startindex);

% Check if we should recurse
if length(qdr_con)>2*maxnsocp
    % recursing costs, so do 4 in one step
    ind = 1+ceil(length(qdr_con)*(0:0.25:1));
    [F_struc1,K1,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(ind(1):ind(2)-1),nvars,maxnsocp,startindex+ind(1)-1);
    [F_struc2,K2,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(ind(2):ind(3)-1),nvars,maxnsocp,startindex+ind(2)-1);
    [F_struc3,K3,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(ind(3):ind(4)-1),nvars,maxnsocp,startindex+ind(3)-1);
    [F_struc4,K4,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(ind(4):ind(5)-1),nvars,maxnsocp,startindex+ind(4)-1);
    F_struc = [F_struc F_struc1 F_struc2 F_struc3 F_struc4];
    K.q = [K1.q K2.q K3.q K4.q];
    K.q(K.q==0)=[];
    return
elseif length(qdr_con)>maxnsocp
    mid = ceil(length(qdr_con)/2);
    [F_struc1,K1,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(1:mid),nvars,maxnsocp,startindex);
    [F_struc2,K2,KCut] = recursive_socp_fix(F,[],K,KCut,qdr_con(mid+1:end),nvars,maxnsocp,startindex+mid);
    F_struc = [F_struc  F_struc1  F_struc2];
    K.q = [K1.q K2.q];
    K.q(K.q==0)=[];
    return
end

% second order cone constraints
for i = 1:length(qdr_con)
    constraints = qdr_con(i);
    
    [n,m] = size(F.clauses{constraints}.data);
    ntimesm = n*m; %Just as well pre-calc
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    
    data = getbase(F.clauses{constraints}.data);
    if isreal(data)
        mapX = [1 1+lmi_variables];
        [ix,jx,sx] = find(data);
        F_structemp = sparse(mapX(jx),ix,sx,1+nvars,ntimesm);
    else
        n = n+(n-1);
        ntimesm = n*m;
        F_structemp  = spalloc(ntimesm,1+nvars,0);
        data = [data(1,:);real(data(2:end,:));imag(data(2:end,:))];
        F_structemp(:,[1 1+lmi_variables(:)'])= data;
        F_structemp = F_structemp';
    end
    % ...and add them together (efficient for large structures)
    F_struc = [F_struc F_structemp];
   % K.q(i+startindex-1) = n;
    K.q = [K.q ones(1,m)*n];
end
K.q(K.q==0)=[];

function [F_struc,K,KCut] = recursive_msocp_fix(F,F_struc,K,KCut,qdr_con,nvars,maxnsocp,startindex);

if isequal(K.q,0)
    K.q = [];
end
% second order cone constraints
for i = 1:length(qdr_con)
    constraints = qdr_con(i);
    
    [n,m] = size(F.clauses{constraints}.data);
    ntimesm = n*m; %Just as well pre-calc
    
    % Which variables are needed in this constraint
    lmi_variables = getvariables(F.clauses{constraints}.data);
    
    data = getbase(F.clauses{constraints}.data);
    if isreal(data)
        mapX = [1 1+lmi_variables];
        [ix,jx,sx] = find(data);
        F_structemp = sparse(mapX(jx),ix,sx,1+nvars,ntimesm);
    else
        n = n+(n-1);
        ntimesm = n*m;
        F_structemp  = spalloc(ntimesm,1+nvars,0);
        data = [data(1,:);real(data(2:end,:));imag(data(2:end,:))];
        F_structemp(:,[1 1+lmi_variables(:)'])= data;
        F_structemp = F_structemp';
    end
    % ...and add them together (efficient for large structures)
    F_struc = [F_struc F_structemp];
    K.q = [K.q ones(1,m)*n];
end






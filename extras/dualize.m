function [Fdual,objdual,X,t,err,complexInfo] = dualize(F,obj,auto,extlp,extend,options)
% DUALIZE Create the dual of an SDP given in primal form
%
% [Fd,objd,X,t,err] = dualize(F,obj,auto)
%
% Input
%  F    : Primal constraint in form AX=b+dt, X>0, t free.
%  obj  : Primal cost CX+ct
%  auto : If set to 0, YALMIP will not automatically handle variables
%        and update variable values when the dual problem is solved.
%  extlp: If set to 0, YALMIP will not try to perform variables changes in
%         order to convert simple translated LP cones (as in x>1) to
%         standard unit cone constraints (x>0)
%
% Output
%  Fd  : Dual constraints in form C-A'y>0, c-dy==0
%  obj : Dual cost b'y (to be MAXIMIZED!)
%  X   : The detected primal cone variables
%  t   : The detected primal free variables
%  err : Error status (returns 0 if no problems)
%
% Example
%  See the HTML help.
%
% See also DUAL, SOLVESDP, PRIMALIZE

% Check for unsupported problems

if isempty(F)
    F = ([]);
end

complexInfo = [];

LogDetTerm = 0;
if nargin < 2
    obj = [];
elseif isa(obj,'logdet')
    Plogdet  = getP(obj);
    gainlogdet = getgain(obj);
    if ~all(gainlogdet <= 0)
        error('There are nonconvex logdet terms in the objective')
    end
    if ~all(gainlogdet == -1)
        error('DUALIZE does currently not support coefficients before LOGDET terms')
    end
    obj = getcx(obj);
    if ~isempty(obj)
        if ~is(obj,'linear')
            error('DUALIZE does not support nonlinear terms in objective (except logdet terms)')
        end
    end
    LogDetTerm = 1;
    Ftemp = ([]);
    for i = 1:length(Plogdet)
        Ftemp = Ftemp + [(Plogdet{i} >= 0) : 'LOGDET'];
    end
    F = Ftemp + F;
end

err = 0;
p1 = ~isreal(obj);%~(isreal(F) & isreal(obj));
p2 = ~(islinear(F) & islinear(obj));
p3 = any(is(F,'integer')) | any(is(F,'binary'));
if p1 | p2 | p3
    if nargout == 5
        Fdual = ([]);objdual = [];y = []; X = []; t = []; err = 1;
    else
        problems = {'Cannot dualize complex-valued problems','Cannot dualize nonlinear problems','Cannot dualize discrete problems'};
        error(problems{min(find([p1 p2 p3]))});
    end
end

if nargin<5 || isempty(extend)
    extend = 1;
end

if extend    
    if nargin < 6 || isempty(options)    
        options = sdpsettings;
    end
    options.dualize = 1;
    options.allowmilp = 0;
    options.solver = '';
    [F,failure,cause] = expandmodel(F,obj,options);
    if failure
        error('Failed during convexity propagation. Avoid nonlinear operators when applying dualization.');
    end
end

if nargin<3 || isempty(auto)
    auto = 1;
end

if nargin<4 || isempty(extlp)
    extlp = 1;
end

% Cones and equalities
F_AXb = ([]);


% Shiftmatrix is a bit messy at the moment.
% We want to be able to allow cones X>shift
% by using a new variable X-shift = Z
shiftMatrix = {};
X={};

% First, get variables in initial SDP cones
% We need this to avoid adding the same variable twice
% when we add simple LP constraints (as in P>0, P(1,3)>0)
varSDP = [];
SDPset = zeros(length(F),1);
ComplexSDPset = zeros(length(F),1);
isSDP = is(F,'sdp');
for i = 1:length(F)
    if isSDP(i);        
        Fi = lmi2sdpvar(F,i);
        if is(Fi,'shiftsdpcone')
            vars = getvariables(Fi);
            if isempty(findrows(varSDP,[vars(1) vars(end)]))
                SDPset(i) = 1;
                varSDP = [varSDP;vars(1) vars(end)];
                shiftMatrix{end+1} = getbasematrix(Fi,0);
                X{end+1}=Fi;
                if is(Fi,'complex')
                    ComplexSDPset(i) = 1;
                end
            end
        end
    end
end
F_SDP = F(find(SDPset));
F = F(find(~SDPset));

% Same thing for second order cones
% However, we must not add any SOC cones
% that we already defined as SDP cones
varSOC = [];
SOCset = zeros(length(F),1);
isSOCP = is(F,'socp');
for i = 1:length(F)
    if isSOCP(i);%is(F(i),'socp')
        Fi = sdpvar(F(i));
        if is(Fi,'socone')
            vars = getvariables(Fi);
            % Make sure these variables are not SDP cone variables
            % This can actually only happen for (X>0) + (Xcone((2:end,1),X(1)))
            if ~isempty(varSDP)
                inSDP = any(varSDP(:,1)<=vars(1)& vars(1) <=varSDP(:,2)) | any(varSDP(:,1)<=vars(end)& vars(end) <=varSDP(:,2));
            else
                inSDP = 0;
            end
            if ~inSDP
                SOCset(i) = 1;
                vars = getvariables(Fi);
                varSOC = [varSOC;vars(1) vars(end)];
                shiftMatrix{end+1} = getbasematrix(Fi,0);
                X{end+1}=Fi;
            end
        end
    end
end
F_SOC = F(find(SOCset));
F = F(find(~SOCset));

% Merge SDP and SOC data
varCONE = [varSDP;varSOC];
F_CONE = F_SDP + F_SOC;

% Detect primal-variables from SDP diagonals (not tested substantially...)
implicit_positive = detect_diagonal_terms(F);

% Find all LP constraints, add slacks and extract simple cones
% to speed up things, we treat LP cone somewhat different
% compared to the conceptually similiar SOCP/SDP cones
% This code is pretty messy, since there are a lot off odd
% cases to take care of (x>0 and x>1 etc etc)
elementwise = is(F,'element-wise');
elementwise_index = find(elementwise);
if ~isempty(elementwise_index)

    % Find element-wise inequalities
    Flp = F(elementwise_index);
    % Add constraint originating from diagonals in dual-type LMIs
    if ~isempty(implicit_positive)        
        implicit_positive = setdiff(implicit_positive,getvariables(F_CONE));
        if ~isempty(implicit_positive) 
            Flp = Flp + (recover(implicit_positive) >= 0);
        end
    end

    F = F(find(~elementwise)); % remove these LPs
    % Find LP cones
    lpconstraint = [];
    for i = 1:length(Flp)
        temp = sdpvar(Flp(i));
        if min(size(temp))>1
            temp = temp(:);
        end
        lpconstraint = [lpconstraint reshape(temp,1,length(temp))];
    end

    % Find all constraints of type a_i+x_i >0 and extract the unique and
    % most constraining inequalities (i.e. remove redundant lower bounds)
    base = getbase(lpconstraint);
    temp = sum(base(:,2:end)~=0,2)==1;
    candidates = find(temp);
    if length(candidates)>0
        % The other ones...
        alwayskeep = find(sum(base(:,2:end)~=0,2)~=1);
        w1 =  lpconstraint(alwayskeep);
        if all(temp)
            w2 = lpconstraint;
        else
            w2 =  lpconstraint(candidates);
        end

        % Find unique rows
        base = getbase(w2);
        [i,uniquerows,k] = unique(base(:,2:end)*randn(size(base,2)-1,1));
        aUniqueRow=k(:)';
        keep = [];
        rhsLP = base(:,1);
        rr = histc(k,[1:max(k)]);
        if all(rr==1)
            lpconstraint = [w1 w2];
        else
            for i=1:length(k)
                sameRow=find(k==k(i));
                if length(sameRow)==1
                    keep = [keep sameRow];
                else
                    rhs=base(sameRow,1);
                    [val,pos]=min(rhsLP(sameRow));
                    keep = [keep sameRow(pos)];
                end
            end
            lpconstraint = [w1 w2(unique(keep))];
        end

    end

    % LP cone will be saved in a vector for speed
    x = [];

    % Pure cones of the type x>0
    base = getbase(lpconstraint);
    purelpcones = (base(:,1)==0) & (sum(abs(base(:,2:end)),2)==1) & (sum(base(:,2:end)==1,2)==1);
    if ~isempty(purelpcones)
        if all(purelpcones)
            x  = [x lpconstraint];
        else
            x  = [x lpconstraint(purelpcones)];
        end
        lpconstraint = lpconstraint(find(~purelpcones));
    end


    % Translated cones x>k, k positive
    % User does not want to make variable changes based on k
    % But if k>=0, we can at-least mark x as a simple LP cone variable and
    % thus avoid a free variable.
    if ~extlp & ~isempty(lpconstraint)
        base = getbase(lpconstraint);
        lpcones = (base(:,1)<0) & (sum(abs(base(:,2:end)),2)==1) & (sum(base(:,2:end)==1,2)==1);
        if ~isempty(find(lpcones))
            s = recover(getvariables(lpconstraint(find(lpcones))));
            x  = [x reshape(s,1,length(s))];
        end
    end

    % Translated cones x>k
    % Extract these and perform the associated variable change y=x-k
    if ~isempty(lpconstraint)%Avoid warning in 5.3.1
        base = getbase(lpconstraint);
        lpcones = (sum(abs(base(:,2:end)),2)==1) & (sum(base(:,2:end)==1,2)==1);
        if ~isempty(lpcones) & extlp
            x  = [x lpconstraint(find(lpcones))];
            nlp = lpconstraint(find(~lpcones));
            if ~isempty(nlp)
                s = sdpvar(1,length(nlp));
                F_AXb = F_AXb + (nlp-s==0);
                x = [x s];
            end
        elseif length(lpconstraint) > 0
            s = sdpvar(1,length(lpconstraint));
            x = [x s]; % New LP cones
            F_AXb = F_AXb + (lpconstraint-s==0);
        end
    end


    % Sort asccording to variable index
    % (Code below assumes x is sorted in increasing variables indicies)
    base = getbase(x);base = base(:,2:end);[i,j,k] = find(base);
    if ~isequal(i,(1:length(x))')
        x = x(i);
    end
    xv = getvariables(x);

    % For mixed LP/SDP problems, we must ensure that LP cone variables are
    % not actually an element in a SDP cone variable
    if ~isempty(varCONE)
        keep = zeros(length(x),1);
        for i = 1:length(xv)
            if any(varCONE(:,1)<=xv(i) & xv(i) <=varCONE(:,2))
            else
                keep(i) = 1;
            end
        end
        if ~all(keep)
            % We need to add some explicit constraints on some elements and
            % remove the x variables since they are already in a cone
            % variable
            xcone = x(find(~keep));
            s = sdpvar(1,length(xcone));
            F_AXb = F_AXb + (xcone-s==0);
            x = x(find(keep));
            x = [x s];
        end
    end
else
    x = [];
end

% Check for mixed cones, ie terms C-A'y > 0.
keep = ones(length(F),1);
isSDP  = is(F,'sdp');
isSOCP = is(F,'socp');
isVecSOCP = is(F,'vecsocp');

% Pre-allocate all SDP slacks in one call
% This is a lot faster
if nnz(isSDP) > 0
    SDPindicies = find(isSDP)';
    for i = 1:length(SDPindicies)%find(isSDP)'
        Fi = lmi2sdpvar(F,SDPindicies(i));
       % Fi = sdpvar(F(SDPindicies(i)));
        ns(i) = size(Fi,1);
        ms(i) = ns(i);
        isc(i) = is(Fi,'complex');
    end
    if any(isc)
        for i = 1:length(ns)
            if isc(i)
                Slacks{i} = sdpvar(ns(i),ns(i),'hermitian','complex');
            else
                Slacks{i} = sdpvar(ns(i),ns(i));
            end
        end
    else
        Slacks = sdpvar(ns,ms);
    end
    if ~isa(Slacks,'cell')
        Slacks = {Slacks};
    end
end
prei = 1;
for i = 1:length(F)
    if isSDP(i)
        % Semidefinite dual-form cone
        %Fi = sdpvar(F(i));
        % MUch faster low-level special function
        Fi = lmi2sdpvar(F,i);
        n  = size(Fi,1);
        %      S  = sdpvar(n,n);
        S  = Slacks{prei};prei = prei + 1;
        slack = Fi-S;
        ind = find(triu(reshape(1:n^2,n,n)));
        if is(slack,'complex')
            F_AXb =  F_AXb + (real(slack(ind))==0) + (imag(slack(ind))==0);
        else
            F_AXb =  F_AXb + (slack(ind)==0);
        end
        F_CONE = F_CONE + lmi(S,[],[],[],1);
        shiftMatrix{end+1} = spalloc(n,n,0);
        X{end+1}=S;
        keep(i)=0;
    elseif isSOCP(i)
        % SOC dual-form cone
        Fi = sdpvar(F(i));
        n  = size(Fi,1);
        S  = sdpvar(n,1);
        %        S = Slacks{i};
        slack = Fi-S;
        if is(slack,'complex')
            F_AXb =  F_AXb + (real(slack)==0) + (imag(slack)==0);
        else
            F_AXb =  F_AXb + (slack==0);
        end
        F_CONE = F_CONE + (cone(S(2:end),S(1)));
        shiftMatrix{end+1} = spalloc(n,1,0);
        X{end+1}=S;
        keep(i)=0;
    elseif isVecSOCP(i)
         % Vectorized SOC dual-form cone
        Fi = sdpvar(F(i));
        [n,m]  = size(Fi);
        S  = sdpvar(n,m,'full');        
        slack = Fi-S;
        if is(slack,'complex')
            F_AXb =  F_AXb + (real(slack)==0) + (imag(slack)==0);
        else
            F_AXb =  F_AXb + (slack==0);
        end
        F_CONE = F_CONE + (cone(S));
        shiftMatrix{end+1} = spalloc(n,m,0);
        X{end+1}=S;
        keep(i)=0;
    end
end

% Now, remove all mixed cones...
F = F(find(keep));

% Get the equalities
AXbset = is(F,'equality');
if any(AXbset)
    % Get the constraints
    F_AXb = F_AXb + F(find(AXbset));
    complex = find(is(F_AXb,'complex'));
    if ~isempty(complex)
        F_AXb_complex = F_AXb(complex);
        F_AXb(complex)=[]; 
        rEQ = real(sdpvar(F_AXb_complex));
        iEQ = imag(sdpvar(F_AXb_complex));
        if ~isempty(rEQ) && isa(rEQ,'sdpvar')
            F_AXb = F_AXb + (rEQ == 0);
        end
        if ~isempty(iEQ) && isa(iEQ,'sdpvar')
            F_AXb = F_AXb + (iEQ == 0);
        end
    end
    F = F(find(~AXbset));
end

% Is there something we missed in our tests?
if length(F)>0
    error('DUALIZE can only treat standard SDPs (and LPs) at the moment.')
end

% If complex SDP cone, we reformulate and call again on a real-valued
% problem. This leads to twice the amount of work, but it is a quick fix
% for the moment
if any(is(F_CONE,'complexsdpcone'))
    F_NEWCONES = [];
    top = 1;
    for i = 1:length(X)      
        if is(X{i},'complexsdpcone')
            Xreplace{top} = X{i};
            n = length(X{i});
            Xnew{top} = sdpvar(2*n);
            
            
            rQ = real(Xreplace{top});
            iQ = imag( Xreplace{top});
            L1 = Xnew{top}(1:n,1:n);
            L3 = Xnew{top}(n+1:end,n+1:end);
            L2 = Xnew{top}(1:n,n + 1:end);
                        
            s0r = getvariables(rQ);
            s1r = getvariables(L1);
            s2r = getvariables(L3);
            r0 = recover(s0r);
            r1 = recover(s1r);
            r2 = recover(s2r);
            
            s0i = getvariables(iQ);
            s1i = getvariables(triu(L2,1))';
            s2i = getvectorvariables(L2(find(tril(ones(length(L2)),-1))));
            i0 = recover(s0i);
            i1 = recover(s1i);
            i2 = recover(s2i);
           
            replacement = [r1+r2;i1-i2];
            if ~isempty(F_AXb)
                F_AXb = remap(F_AXb,[s0r s0i],replacement);
            end
            if ~isempty(obj)
                obj = remap(obj,[s0r s0i],replacement);
            end
                                   
            X{i} = Xnew{top};
            top = top + 1;
        end
        if is(X{i},'hermitian')
            F_NEWCONES = [F_NEWCONES, X{i} >= 0];
        else 
            F_NEWCONES = [F_NEWCONES, cone(X{i})];
        end
    end  
    F_reformulated = [F_NEWCONES, F_AXb, x>=0];
    complexInfo.replaced = Xreplace;
    complexInfo.new = Xnew;
    [Fdual,objdual,X,t,err] = dualize(F_reformulated,obj,auto,extlp,extend);
    return
end

% Sort the SDP cone variables X according to YALMIP
% This is just to simplify some indexing later
ns = [];
first_var = [];
for i = 1:length(F_CONE)
    ns = [ns length(X{i})];
    first_var = [first_var min(getvariables(X{i}))];
end
[sorted,index] = sort(first_var);
X={X{index}};
shiftMatrix={shiftMatrix{index}};

shift = [];
for i = 1:length(F_CONE)
    ns(i) = length(X{i});
    if size(X{i},2)==1 | (size(X{i},1) ~= size(X{i},2))
        % SOCP
        shift = [shift;shiftMatrix{i}(:)];
    else
        % SDP
        ind =  find(tril(reshape(1:ns(i)^2,ns(i),ns(i))));
        shift = [shift;shiftMatrix{i}(ind)];
    end
end

% Free variables (here called t) is everything except the cone variables
X_variables = getvariables(F_CONE);
x_variables = getvariables(x);
Xx_variables = [X_variables x_variables];

other_variables = [getvariables(obj) getvariables(F_AXb)];
% For quadratic case
%other_variables = [depends(obj) getvariables(F_AXb)];

all_variables = uniquestripped([other_variables Xx_variables]);

% Avoid set-diff
if isequal(all_variables,Xx_variables)
    t_variables = [];
    t_in_all = [];
    t = [];
else
    t_variables = setdiff(all_variables,Xx_variables);
    ind = ismembcYALMIP(all_variables,t_variables);
    t_in_all = find(ind);
    t = recover(t_variables);
end

ind = ismembcYALMIP(all_variables,x_variables);
x_in_all = find(ind);
ind = ismembcYALMIP(all_variables,X_variables);
X_in_all = find(ind);

vecF1 = [];
nvars = length(all_variables);
for i = 1:length(F_AXb)
    AXb = sdpvar(F_AXb(i));
    mapper = find(ismembcYALMIP(all_variables,getvariables(AXb)));

    [n,m] = size(AXb);
    data = getbase(AXb);
    [iF,jF,sF] = find(data);
    if 1 % New
        smapper = [1 1+mapper(:)'];
        F_structemp = sparse(iF,smapper(jF),sF,n*m,1+nvars);
    else
        F_structemp  = spalloc(n*m,1+nvars,nnz(data));
        F_structemp(:,[1 1+mapper(:)'])= data;
    end
    vecF1 = [vecF1 F_structemp'];
end
vecF1 = vecF1';

%Remove trivially redundant constraints
h = 1+rand(size(vecF1,2),1);
h = vecF1*h;
% INTVAL possibility
%[dummy,uniquerows] =  uniquesafe(h);
[dummy,uniquerows] =  uniquesafe(mid(h));
if length(uniquerows) < length(h)
    % Sort to ensure run-to-run consistency
    vecF1 = vecF1((sort(uniquerows)),:);
end

if isempty(obj)
    vecF1(end+1,1) = 0;
else
    if is(obj,'linear')
        mapper = find(ismembcYALMIP(all_variables,getvariables(obj)));
        [n,m] = size(obj);
        data = getbase(obj);
        [iF,jF,sF] = find(data);
        if 1
            smapper = [1 1+mapper(:)'];
            F_structemp = sparse(iF,smapper(jF),sF,n*m,1+nvars);
        else
            F_structemp  = spalloc(n*m,1+nvars,nnz(data));
            F_structemp(:,[1 1+mapper(:)'])= data;
        end
        vecF1 = [vecF1;F_structemp];
    else
        % FIX: Generalize to QP duality
        % min c'x+0.5x'Qx, Ax==b, x>=0
        % max b'y-0.5x'Qx, c-A'y+Qx >=0
        [Q,c,xreally,info] = quaddecomp(obj,recover(all_variables))
        mapper = find(ismembcYALMIP(all_variables,getvariables(c'*xreally)));
        [n,m] = size(c'*xreally);
        data = getbase(c'*xreally);
        F_structemp  = spalloc(n*m,1+nvars,nnz(data));
        F_structemp(:,[1 1+mapper(:)'])= data;
        vecF1 = [vecF1;F_structemp]
    end

end

vecF1(end+1,1) = 0;
Fbase = vecF1;

%Fbase = unique(Fbase','rows')';

b = Fbase(1:end-2,1);

Fbase = -Fbase(1:end-1,2:end);
vecA = [];
Fbase_t = Fbase(:,t_in_all);
Fbase_x = Fbase(:,x_in_all);
Fbase_X = Fbase;
%Fbase_X(:,unionstripped(t_in_all,x_in_all)) = [];
if 1
    removethese = unique([t_in_all x_in_all]);
    if length(removethese) > 0.5*size(Fbase_X,2)
        Fbase_X = Fbase_X(:,setdiff(1:size(Fbase_X,2),removethese));
    else
        Fbase_X(:,[t_in_all x_in_all]) = [];
    end
else
    removecols = uniquestripped([t_in_all x_in_all]);
    if ~isempty(removecols)
        [i,j,k] = find(Fbase_X);
        keep = find(~ismember(j,removecols));
        i = i(keep);
        k = k(keep);
        j = j(keep);
        map = find(1:length(unique(j)),j);

    end
end

% Shift due to translated dual cones X = Z+shift
if length(shift > 0)
    b = b + Fbase_X(1:end-1,:)*shift;
end
if length(x)>0
    % Arrgh
    base = getbase(x);
    constant = base(:,1);
    base = base(:,2:end);[i,j,k] = find(base);
    b = b + Fbase_x(1:end-1,:)*constant(i);
end

start = 0;
n_cones = length(ns);
% All LPs in one cone
if length(x)>0
    n_lp = 1;
else
    n_lp = 0;
end
n_free = length(t_variables);

% SDP cones
for j = 1:1:n_cones

    if size(X{j},1)==size(X{j},2)
        % SDP cone
        ind = reshape(1:ns(j)^2,ns(j),ns(j));
        ind = find(tril(ind));

        % Get non-symmetric constraint AiX=b
        Fi = Fbase_X(1:length(b),start+(1:length(ind)))'/2;

        if 1
            [iF,jF,sF] = find(Fi);
            iA = ind(iF);
            jA = jF;
            sA = sF;
            the_col = 1+floor((iA-1)/ns(j));
            the_row = iA-(the_col-1)*ns(j);
            these_must_be_transposed = find(the_row > the_col);
            if ~isempty(these_must_be_transposed)
                new_rowcols = the_col(these_must_be_transposed) + (the_row(these_must_be_transposed)-1)*ns(j);
                iA = [iA;new_rowcols];
                jA = [jA;jA(these_must_be_transposed)];
                sA = [sA;sA(these_must_be_transposed)];
            end
            % Fix diagonal term
            diags = find(diag(1:ns(j)));
            id = find(ismembcYALMIP(iA,diags));
            sA(id) = 2*sA(id);
            Ai = sparse(iA,jA,sA,ns(j)^2,length(b));

        else % Old slooooooow version
            Ai = spalloc(ns(j)^2,length(b),nnz(Fi));
            Ai(ind,:) = Fi;
            % Symmetrize
            [rowcols,varindicies,vals]=find(Ai);
            the_col = 1+floor((rowcols-1)/ns(j));
            the_row = rowcols-(the_col-1)*ns(j);
            these_must_be_transposed = find(the_row > the_col);
            if ~isempty(these_must_be_transposed)
                new_rowcols = the_col(these_must_be_transposed) + (the_row(these_must_be_transposed)-1)*ns(j);
                Ai(sub2ind(size(Ai),new_rowcols,varindicies(these_must_be_transposed))) = vals(these_must_be_transposed);
            end

            % Fix diagonal term
            diags = find(diag(1:ns(j)));
            Ai(diags,:) = 2*Ai(diags,:);
        end

        %        if norm(Ai-Ai2,inf)>1e-12
        %            error
        %        end


        vecA{j} = Ai;

        start = start + length(ind);
    else
        % Second order cone
        ind = 1:prod(size(X{j}));
        %ind = 1:ns(j);

        % Get constraint AiX=b
        Fi = Fbase_X(1:length(b),start+(1:length(ind)))';
        %Ai = spalloc(ns(j),length(b),nnz(Fi));
        Ai = spalloc(prod(size(X{j})),length(b),nnz(Fi));
        Ai(ind,:) = Fi;
        vecA{j} = Ai;
        start = start + length(ind);
    end

end
% LP Cone
if n_lp>0
    Alp=Fbase_x(1:length(b),:)';
end
% FREE VARIABLES
start = 0;
if n_free>0
    Afree=Fbase_t(1:length(b),:)';
end

% COST MATRIX
% SDP CONE
start = 0;
for j = 1:1:n_cones
    if size(X{j},1)==size(X{j},2)
        %ind = reshape(1:ns(j)^2,ns(j),ns(j));
        %ind = find(tril(ind));
        %C{j} = spalloc(ns(j),ns(j),0);
        %C{j}(ind) = -Fbase_X(end,start+(1:length(ind)));
        %C{j} = (C{j}+ C{j}')/2;
        %start = start + length(ind);

        ind = reshape(1:ns(j)^2,ns(j),ns(j));
        [indi,indj] = find(tril(ind));
        C{j} = sparse(indi,indj,-Fbase_X(end,start+(1:length(indi))),ns(j),ns(j));
        C{j} = (C{j}+ C{j}')/2;
        start = start + length(indi);
    else
        %ind = 1:ns(j);
        ind = 1:prod(size(X{j}));
        C{j} = spalloc(ns(j),1,0);
        C{j}(ind) = -Fbase_X(end,start+(1:length(ind)));
        start = start + length(ind);
    end
end
% LP CONE
for j = 1:1:n_lp
    Clp = -Fbase_x(end,:)';
end
% FREE CONE
if n_free>0
    Cfree = -Fbase_t(end,1:end)';
end

% Create dual
if length(b) == 0
    error('Dual program is somehow trivial (0 variables in dual)');
end
y = sdpvar(length(b),1);
yvars = getvariables(y);
cf = [];
Af = [];
Fdual = ([]);
for j = 1:n_cones
    if size(X{j},1)==size(X{j},2)
        % Yep, this is optimized...
        S = sdpvar(ns(j),ns(j),[],yvars,[C{j}(:) -vecA{j}]);
        % Fast call avoids symmetry check
        Fdual = Fdual + lmi(S,[],[],[],1);
    else
        Ay = reshape(vecA{j}*y,[],1);
        %Ay = reshape(vecA{j}*y,ns(j),1);
        S = C{j}-Ay;
        S = reshape(S,size(X{j},1),[]);
        %Fdual = Fdual + lmi(cone(S(2:end),S(1)));
        Fdual = Fdual + lmi(cone(S));
    end
end
if n_lp > 0
    keep = any(Alp,2);
    if ~all(keep)
        % Fix for unused primal cones
        preset=find(~keep);
        xfix = x(preset);
        assign(xfix(:),Clp(preset(:)));
    end
    keep = find(keep);
    if ~isempty(keep)
        z = Clp(keep)-Alp(keep,:)*y;
        if isa(z,'double')
            assign(x(:),z(:));
        else
            Fdual = Fdual + lmi(z);
            if ~isequal(keep,(1:length(x))')
                x = x(keep);
            end
            X{end+1} = x(:); % We have worked with a row for performance reasons
        end
    end
end
if n_free > 0
    CfreeAfreey = Cfree-Afree*y;
    if isa(CfreeAfreey,'double')
        if nnz(CfreeAfreey)>0
            error('Trivially unbounded!');
        end
    else
        Fdual = Fdual + (0 == CfreeAfreey);
    end
end

objdual = b'*y;
if auto
    for i = 1:length(X)
        yalmip('associatedual',getlmiid(Fdual(i)),X{i});
    end
    if n_free>0
        yalmip('associatedual',getlmiid(Fdual(end)),t);
    end
end

if LogDetTerm
    for i = 1:length(Plogdet)
        objdual = objdual + logdet(sdpvar(Fdual(i))) + length(sdpvar(Fdual(1)));
    end
    Fdual = Fdual - Fdual(1:length(Plogdet));
end

Fdual = setdualize(Fdual,1);

function implicit_positive = detect_diagonal_terms(F)
F = F(find(is(F,'sdp')));
implicit_positive = [];
for i = 1:length(F)
    Fi = lmi2sdpvar(F,i);
   % Fi = sdpvar(F(i));
    B = getbase(Fi);
    n = sqrt(size(B,1));
    d = 1:(n+1):n^2;
    B = B(d,:);
    candidates = find((B(:,1) == 0) & (sum(B | B,2) == 1) & (sum(B,2) == 1));
    if ~isempty(candidates)
        vars = getvariables(Fi);
        [ii,jj,kk] = find(B(candidates,2:end)');
        ii = ii(:)';
        implicit_positive = [implicit_positive vars(ii)];
    end
end



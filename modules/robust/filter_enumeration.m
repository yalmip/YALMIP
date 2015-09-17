function [F,mptmissing] = filter_enumeration(F_xw,Zmodel,x,w,ops,uncertaintyTypes,separatedZmodel,VariableType)

mptmissing = 0;
if length(F_xw) == 0
    F = [];
    return;
else
    
    if any(Zmodel.K.q) | any(Zmodel.K.s)
        error('Only polytope uncertainty supported in duality based robustification');
    else
        if isempty(intersect(depends(F_xw),getvariables(w)))
            F = F_xw;
        elseif length(uncertaintyTypes)==1 & isequal(uncertaintyTypes{1},'inf-norm')
            
            if any(isinf((separatedZmodel{1}.lb))) | any(isinf(separatedZmodel{1}.ub))
                error('You have unbounded uncertain variables')
            else
                n = length(separatedZmodel{1}.lb);
                vertices = [];
                lb = separatedZmodel{1}.lb(:)';
                ub = separatedZmodel{1}.ub(:)';
                E = dec2bin(0:2^n-1,n)';                
                E = double(E(:))-48;                
                E = reshape(E,n,2^n);
                vertices = (repmat(lb(:),1,2^n) + E.*(repmat(ub(:),1,2^n)-repmat(lb(:),1,2^n)))';
                %for i = 0:2^n-1
                %    vertices = [vertices;lb+dec2decbin(i,n).*(ub-lb)];
                %end
                if ops.verbose
                    disp([' - Enumerated ' num2str(2^n) ' vertices'])
                end
                vertices = unique(vertices,'rows');
                if ops.verbose & 2^n > size(vertices,1)
                    disp([' - Reduced to ' num2str( size(vertices,1)) ' unique vertices'])
                end
                F = replaceVertices(F_xw,w,vertices',VariableType,ops);
            end
            
        elseif length(uncertaintyTypes)==1 & isequal(uncertaintyTypes{1},'simplex')
            
            k = abs(Zmodel.F_struc(1,1));
            n = length(w);
            
            vertices = zeros(n,1);
            for i = 1:n
                v = zeros(n,1);
                v(i) = k;
                vertices = [vertices v];
            end
            if ops.verbose
                disp([' - Enumerated ' num2str(n) ' vertices'])
            end
            vertices = pruneequalities(vertices,Zmodel);
            F = replaceVertices(F_xw,w,vertices,VariableType,ops);
            
            
        else
            % FIX : Assumes all uncertainty in all constraints
            K = Zmodel.K;
            A = -Zmodel.F_struc((1+K.f):(K.f + K.l),2:end);
            b =  Zmodel.F_struc((1+K.f):(K.f + K.l),1);
            try
                % Some preprocessing to extract bounds from equality
                % constraints in order to make the uncertainty polytope
                % bounded (required since we are going to run vertex
                % enumeration)
                % We might have x>=0, sum(x)=1, and this code simply extracts
                % the implied bounds x<=1
                [lo,up] = findulb(Zmodel.F_struc(1:K.f + K.l,:),K);
                Zmodel.lb = lo;Zmodel.ub = up;
                Zmodel = propagate_bounds_from_equalities(Zmodel);
                up = Zmodel.ub;
                lo = Zmodel.lb;
                upfi = find(~isinf(up));
                lofi = find(~isinf(lo));
                aux = Zmodel;
                aux.F_struc = [aux.F_struc;-lo(lofi) sparse(1:length(lofi),lofi,1,length(lofi),size(A,2))];
                aux.F_struc = [aux.F_struc;up(upfi) -sparse(1:length(upfi),upfi,1,length(upfi),size(A,2))] ;
                aux.K.l = aux.K.l + length(lofi) + length(upfi);
                K = aux.K;
                A = -aux.F_struc((1+K.f):(K.f + K.l),2:end);
                b =  aux.F_struc((1+K.f):(K.f + K.l),1);
                P = polytope(full(A),full(b));
                try
                    vertices = extreme(P)';
                catch
                    error('The uncertainty space is unbounded (could be an artefact of YALMIPs modelling of nonolinear oeprators).')
                end
                %if ~isbounded(P)
                %    error('The uncertainty space is unbounded (could be an artefact of YALMIPs modelling of nonolinear oeprators).')
                %else
                %    vertices = extreme(polytope(A,b))';
                %end
            catch
                mptmissing = 1;
                if ops.verbose>0
                    %lasterr
                    disp(' - Enumeration of uncertainty polytope failed. Missing Multiparametric Toolbox?')
                    disp(' - Switching to duality based approach')
                    %disp('You probably need to install MPT (needed for vertex enumeration)')
                    %disp('http://control.ee.ethz.ch/~joloef/wiki/pmwiki.php?n=Solvers.MPT')
                    %disp('Alternatively, you need to add bounds on your uncertainty.')
                    %disp('Trying to switch to dualization approach')
                    %error('MPT missing');
                end
                F = [];
                return
            end
          
            % The vertex enumeration was done without any equality constraints.
            % We know check all vertices so see if they satisfy equalities.
            vertices = pruneequalities(vertices,Zmodel);
            if ops.verbose
                disp([' - Enumerated ' num2str(size(vertices,2)) ' vertices'])
            end           
            F = replaceVertices(F_xw,w,vertices,VariableType,ops);            
        end
    end
end


function F = replaceVertices(F_xw,w,vertices,VariableType,ops)
% Doing LP constraints in a vectorized manner saves a lot of time

F_xw_lp = F_xw(find(is(F_xw,'elementwise')));
F_xw_socp_sdp = F_xw -  F_xw_lp;
F = ([]);

x_Flp = depends(F_xw_lp);
uncAux = yalmip('auxvariablesW');
uncAux = recover(intersect(x_Flp,VariableType.aux_with_w_dependence));
if isequal(ops.robust.auxreduce,'none')
    uncAux = [];
end
    
w = flush(w);

if length(F_xw_lp)>0
    rLP = [];
    if ~isempty(uncAux)     
        z = sdpvar(repmat(length(uncAux),1,size(vertices,2)),repmat(1,1,size(vertices,2)),'full');
    end
    for i = 1:size(vertices,2)
        temp = replace(sdpvar(F_xw_lp),w,vertices(:,i),0);
        if ~isempty(uncAux)            
            temp = replace(temp,uncAux,z{i});
        end
        rLP = [rLP;temp];
    end
    
    % FIXME: More general detection of silly constraints
    if isa(rLP,'double') & all(rLP>=-eps^0.75)
        F = ([]);
    else
        % Easily generates redundant constraints
        [aux,index] = uniquesafe(getbase(rLP),'rows');
        try
            F = (rLP(index(randperm(length(index)))) >= 0);
        catch
            1
        end
    end
end

% Remaining conic stuff
for j = 1:length(F_xw_socp_sdp)
    for i = 1:size(vertices,2)
        temp = replace(F_xw_socp_sdp(j),w,vertices(:,i),0);
        if ~isempty(uncAux)
            temp = replace(temp,uncAux,z{i});
        end
        F = F + lmi(temp);
    end
end

function vertices = pruneequalities(vertices,Zmodel)
K = Zmodel.K;
% The vertex enumeration was done without any equality constraints.
% We know check all vertices so see if they satisfy equalities.
if K.f > 0
    Aeq = -Zmodel.F_struc(1:K.f,2:end);
    beq =  Zmodel.F_struc(1:K.f,1);
    feasible = sum(abs(Aeq*vertices - repmat(beq,1,size(vertices,2))),1) < 1e-6;
    vertices = vertices(:,feasible);
    if isempty(feasible)
        error('The uncertainty space is infeasible.')
    end
end

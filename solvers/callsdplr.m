function output = callsdplr(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
lowrankdetails =  interfacedata.lowrankdetails;

% Create the parameter structure
pars = options.sdplr;
pars.printlevel = options.verbose;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% rmfield is slow...
Knew.l = K.l;
Knew.s = K.s;
K = Knew;

if options.savedebug
    save sdplrdebug F_struc c K pars -V6
end

% *********************************************
%FIND LOW RANK STRUCTURES
% *********************************************
problem = 0;
lrA = [];
if ~isempty(lowrankdetails) | (options.sdplr.maxrank>0 & (K.s(1)>0))
    showprogress('Detecting low rank data',options.showprogress);
    sdploc = K.l+cumsum([1 K.s.^2]);
    k = 1;
    % FIX : Lazy code copy...
    [ix,jx,sx] = find(F_struc);
    if isempty(lowrankdetails)
        % Okay, this means we have to go through everything...
        % Get all data in F_struc for later
        for lmiid = 1:length(K.s)
            removethese = zeros(1,size(F_struc,2)-1);
            for i = 1:size(F_struc,2)-1
                Fi = reshape(F_struc(sdploc(lmiid):sdploc(lmiid+1)-1,i+1),K.s(lmiid),K.s(lmiid));
                if nnz(Fi)>0
                    [D,V] = getfactors(Fi);
                    if length(D) <= options.sdplr.maxrank
                        lrA(k).cons = i;
                        lrA(k).start = sdploc(lmiid);
                        lrA(k).D = D;
                        lrA(k).V = V;
                        k = k+1;
                        removethese(i) = 1;
                    end
                end
            end
            removethese = find(removethese);
            if ~isempty(removethese)
                these = find((sdploc(lmiid+1)-1>= ix) & (ix>=sdploc(lmiid)) & ismember(jx,1+removethese));
                sx(these) = 0;
            end
        end
    else
        % Just check those constraints declared low-rank by user
        for lrdef = 1:length(lowrankdetails)
            for lrconstraint = 1:length(lowrankdetails{lrdef}.id)
                lmiid = lowrankdetails{lrdef}.id(lrconstraint);
                removethese = zeros(1,size(F_struc,2)-1);
                checkthese = lowrankdetails{lrdef}.variables;
                if isempty(checkthese)
                    checkthese = 1:size(F_struc,2)-1;
                end
                for i = checkthese
                    Fi = reshape(F_struc(sdploc(lmiid):sdploc(lmiid+1)-1,i+1),K.s(lmiid),K.s(lmiid));
                    if nnz(Fi)>0
                        [D,V] = getfactors(Fi);
                        if (options.sdplr.maxrank == 0) | (options.sdplr.maxrank ~= 0 & (length(D) <= options.sdplr.maxrank))
                            lrA(k).cons = i;
                            lrA(k).start = sdploc(lmiid);
                            lrA(k).D = D;
                            lrA(k).V = V;
                            k = k+1;
                            removethese(i) = 1;
                        end
                    end
                end
                removethese = find(removethese);
                if ~isempty(removethese)
                    these = find((sdploc(lmiid+1)-1>= ix) & (ix>=sdploc(lmiid)) & ismember(jx,1+removethese));
                    sx(these) = 0;
                    %F_struc(sdploc(lmiid):sdploc(lmiid+1)-1,1+removethese) = 0;
                end
            end
        end
        F_struc = sparse(ix,jx,sx,size(F_struc,1),size(F_struc,2));
    end
end


% *********************************************
% CALL SDPLR
% *********************************************
if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if isempty(lrA)    
   [x_s,y_s,info] = sdplr(F_struc(:,2:end),c,F_struc(:,1),K);
else   
   % pars.reduce = 0;
    [x_s,y_s,info] = sdplr(F_struc(:,2:end),full(c),full(F_struc(:,1)),K,pars,lrA);
end
solvertime = toc(solvertime);

% YALMIP format
D_struc = x_s;
x = -y_s;

% No error codes currently...
problem = 0;

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.A = F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.b = c;
    solverinput.K = K;
    solverinput.pars = pars;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);

function [D,V] = getfactors(Fi)
if nnz(Fi)>0
    [v,d] = eig(full(Fi));
    d = diag(d);
    keep = find(abs(d)>1e-6);
    V = v(:,keep);
    D = d(keep);
    % lrA(k).cons = i;
    % lrA(k).start = sdploc(j);
    % lrA(k).D = D;
    % lrA(k).V = V;
    % k = k+1;
else
    D = [];
    V = [];
end



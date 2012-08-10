function output = callmpt3(interfacedata)

% Author Johan Löfberg

% Speeds up solving LPs in mpmilp
global MPTOPTIONS
if ~isstruct(MPTOPTIONS)
    mpt_error
end

% Convert
% interfacedata = pwa_linearize(interfacedata);
Matrices = yalmip2mpt(interfacedata);

% Get some MPT options
options = interfacedata.options;
options.mpt.lpsolver = MPTOPTIONS.lpsolver;
options.mpt.milpsolver = MPTOPTIONS.milpsolver;
options.mpt.verbose = options.verbose;

if options.savedebug
    save mptdebug Matrices
end

if options.mp.unbounded
    Matrices = removeExplorationConstraints(Matrices);
end

[dummy,un] = unique([Matrices.G Matrices.E Matrices.W],'rows');
Matrices.G = Matrices.G(un,:);
Matrices.E = Matrices.E(un,:);
Matrices.W = Matrices.W(un,:);

if isempty(Matrices.binary_var_index)

    showprogress('Calling MPT',options.showprogress);
    solvertime = clock;
    if options.mp.presolve
        [Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);
    end        
    
    if any(Matrices.lb(end-Matrices.nx+1:end) == Matrices.ub(end-Matrices.nx+1:end))
        model = [];
    else        
        model = mpt_solvenode(Matrices,Matrices.lb,Matrices.ub,Matrices,[],options);
    end
    solvertime = etime(clock,solvertime);

else  
    % Pre-solve required on binary problems
    options.mp.presolve = 1;

    solvertime = clock;  
    
 %   if Matrices.qp &  options.mp.algorithm == 3
 %        options.mp.algorithm = 1;
 %   end
     
    switch options.mp.algorithm
        case 1
            showprogress('Calling MPT via enumeration',options.showprogress);
            model = mpt_enumeration_mpmilp(Matrices,options);
        case 2
            % Still experimental and just for fun. Not working!
            showprogress('Calling MPT via parametric B&B',options.showprogress);
            model = mpt_parbb(Matrices,options);         
            
       case 3
            showprogress('Calling MPT via delayed enumeration',options.showprogress);
            %Matrices = initialize_binary_equalities(Matrices)           
            [Matrices.SOS,Matrices.SOSVariables] =  mpt_detect_sos(Matrices);
            [Matrices.lb,Matrices.ub] = mpt_detect_and_improve_bounds(Matrices,Matrices.lb,Matrices.ub,Matrices.binary_var_index,options);                                
            model = mpt_de_mpmilp(Matrices,options,[]);            
             
            
            
        otherwise
    end
    solvertime = etime(clock,solvertime);
end

if isempty(model)
    model = {model};
end

if options.verbose
    if ~isempty(model{1})
        if length(model) == 1
            disp(['-> Generated 1 partition.'])            
        else
            disp(['-> Generated ' num2str(length(model)) ' partitions.'])
        end
    end
end

problem = 0;
infostr = yalmiperror(problem,'MPT');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.Matrices = Matrices;
    solverinput.options  = [];
else
    solverinput = [];
end

% Save all data from the solver?
% This always done
if options.savesolveroutput
    solveroutput.model = model;
    solveroutput.U = interfacedata.used_variables(Matrices.free_var);%(Matrices.free_var <= length( interfacedata.used_variables)));
    solveroutput.x = interfacedata.used_variables(Matrices.param_var);
else
    solveroutput = [];
end

% Standard interface
output.Primal      = nan*ones(length(interfacedata.c),1);
output.Dual        = [];
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;

% 
% function M = pwa_linearize(M)
% 
% if any(M.variabletype ~= 0)
%     [lb,ub] = findulb(M.F_struc,M.K);
%     M.lb = lb;
%     M.ub = ub;
%     nonlinear = find(M.variabletype ~= 0);
%     ok = 1;
%     for i = 1:length(nonlinear)
%         if nnz(M.monomtable(nonlinear(i),:)) == 1 
%             j = find(M.monomtable(nonlinear(i),:));
%             if isinf(lb(j)) | isinf(ub(j))
%                 ok = 0;
%             end
%         else
%             ok = 0;
%         end
%     end
%     if ok
%         N = 5;
%         for i = 1:length(nonlinear)
%            j = find(M.monomtable(nonlinear(i),:));
%            x = linspace(lb(j),ub(j),N);
%            y = x.^M.monomtable(nonlinear(i),j);
%            k = [];
%            m = [];
%            for r = 1:N-1
%              km = polyfit(x(r:r+1),y(r:r+1),1);
%              k = [k km(1)];
%              m = [m km(2)];
%            end
%            % Redefine as linear
%            M.monomtable(nonlinear(i),:) = 0;
%            M.monomtable(nonlinear(i),nonlinear(i)) = 1;
%            M.variabletype(nonlinear(i)) = 0;
%            % Add new binaries
%            M.monomtable = blkdiag(M.monomtable,eye(N-1));
%            M.variabletype = [M.variabletype zeros(1,N-1)];           
%            binaries = [length(M.monomtable)-N+2:length(M.monomtable)];           
%            M.binary_variables = [M.binary_variables binaries];
%          %  M.free_var = [M.free_var binaries];
%            M.ub(nonlinear(i)) = max(y);
%            M.lb(nonlinear(i)) = min(y);
%            M.lb(binaries) = 0;
%            M.ub(binaries) = 1;
%                           
%            M.F_struc = [zeros(1,size(M.F_struc,2));M.F_struc];
%            M.F_struc(1,1) = 1;
%            M.F_struc(1,1+binaries) = -1;
%            M.K.f = M.K.f + 1;
%            A = [];
%            b = [];
%            bigM = 2*(max(y)-min(y));
%            for r = 1:N-1
%                A(end+1,binaries(r)) = -bigM;
%                A(end,nonlinear(i)) = -1;
%                A(end,j) = k(r);
%                b(end+1) = m(r)+bigM;
%            end
%            for r = 1:N-1
%                 A(end+1,binaries(r)) = -bigM;
%                 A(end,nonlinear(i)) = 1;
%                 A(end,j) = -k(r);
%                 b(end+1) = bigM-m(r);
%            end
%            bigM = 20;
%            for r = 1:N-1
%                A(end+1,binaries(r)) = -bigM;
%                A(end,j) = -1;
%                b(end+1) = bigM+x(r+1);
%            end
%            for r = 1:N-1
%                A(end+1,binaries(r)) = -bigM;
%                A(end,j) = 1;
%                b(end+1) = bigM-x(r);
%            end
%                                 
%             M.F_struc = [M.F_struc(1:M.K.f,:);b(:) A;M.F_struc(M.K.f+1:end,:)];
%             M.K.l = M.K.l + 1;                     
%         end
%         M.c(length(M.variabletype)) = 0;
%         M.Q(length(M.variabletype),length(M.variabletype)) = 0;
%     end   
% end
% 
% 
% 

function Matrices = initialize_binary_equalities(Matrices)
binary_var_index = Matrices.binary_var_index;
notbinary_var_index = setdiff(1:Matrices.nu,binary_var_index);
% Detect and extract pure binary equalities. Used for simple pruning
nbin = length(binary_var_index);
only_binary = ~any(Matrices.Aeq(:,notbinary_var_index),2);
Matrices.Aeq_bin = Matrices.Aeq(find(only_binary),binary_var_index);
Matrices.beq_bin = Matrices.beq(find(only_binary),:);

function Matrices = removeExplorationConstraints(Matrices);
candidates = find((~any(Matrices.G,2)) & (sum(Matrices.E | Matrices.E,2) == 1));
if ~isempty(candidates)
    Matrices.bndA = -Matrices.E(candidates,:);
    Matrices.bndb = Matrices.W(candidates,:);
    Matrices.G(candidates,:) = [];
    Matrices.E(candidates,:) = [];
    Matrices.W(candidates,:) = [];
end

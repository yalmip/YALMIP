function [UncertainModel,Uncertainty,VariableType,ops,failure] = decomposeUncertain(F,h,w,ops)

failure = 0;
% Do we have any uncertainty declarations variables?
[F,w] = extractUncertain(F,w);
if isempty(w)
    error('There is no uncertainty in the model.');   
end

% Partition the model into
% F_x  : Constraints in decision variables only
% F_w  : The uncertainty description
% F_xw : The uncertain constraints
% Note that this analysis might also declare som of the auxiliary variables
% as simple uncertain variables. It might also create a new objective
% function in order to have all uncertainty in the constraints
F_original = F;
[VariableType,F_x,F_w,F_xw,h] = robust_classify_variables_newest(F,h,ops,w);

if length(F_w)==0
    error('There is no uncertainty description in the model.');
end

if ops.verbose
    dispsilent(ops.verbose,'***** Starting YALMIP robustification module. *********************');
   if length(w)<length(VariableType.w_variables)
       dispsilent(ops.verbose,[' - Detected ' num2str(length(VariableType.w_variables)) ' uncertain variables (' num2str(length(VariableType.w_variables)-length(w)) ' artificial)']);
   else
       dispsilent(ops.verbose,[' - Detected ' num2str(length(w)) ' uncertain variables']);
   end
end

% Integer variables are OK in x, but not in the uncertainty
integervars = [yalmip('binvariables') yalmip('intvariables')];
ind = find(is(F_original,'integer') | is(F_original,'binary'));
if ~isempty(ind)
    integervars = [integervars getvariables(F(ind))];
    if any(ismember(VariableType.w_variables,integervars))
        failure = 1;
        return
    end
end

% Convert quadratic constraints in uncertainty model to SOCPs. This will
% enable us to use duality based removal of uncertainties in linear
% inequalities. We keep information about quadratic expression though, in
% order to use them if possible in 2-norm explicit maximization
F_w = convertquadratics(F_w);

% Convert quadratic constraints in uncertain model to SOCPs. This will
% enable us to use conic/conic robustification
F_xw = convertquadratics(F_xw);

% Export uncertainty model to numerical format
ops.solver = '';
ops.removeequalities = 0;
[aux1,aux2,aux3,Zmodel] = export(F_w-F_w(find(is(F_w,'uncertain'))),[],ops,[],[],1);

if ~isempty(Zmodel)
    if length(Zmodel.c) ~= length(VariableType.w_variables)
        error('Some uncertain variables are unconstrained.')
    end
else
    error('Failed when exporting a model of the uncertainty.')
end

% The uncertainty model is in the full w-space. However, it might
% happen that the uncertainty is separable. Find groups of uncertain
% variables that can/should be treated separately.
uncertaintyGroups = findSeparate(Zmodel);

if ops.verbose
    dispsilent(ops.verbose,[' - Detected ' num2str(length(uncertaintyGroups)) ' independent group(s) of uncertain variables']);
end

% Trying to take care of cases such as norm([x+w;x-w]), i.e. epigraphs with
% uncertainty.
%[F_xw,ops] = prepareforAuxiliaryRemoval(VariableType,F_xw,F_w,ops);
x = recover(VariableType.x_variables);
w = recover(VariableType.w_variables);
ops.robust.forced_enumeration = 0;
switch ops.robust.auxreduce
    case {'none','affine','','projection','enumeration'}
    otherwise
        disp(' ');
        dispsilent(ops.verbose,['The flag ''auxreduce'' is wrong. Turning off removal of auxilliary variables']);
        disp(' ');
        ops.robust.auxreduce = 'none';
end

if ~isempty(VariableType.aux_with_w_dependence)
    if ~strcmp('none',ops.robust.auxreduce)
        dispsilent(ops.verbose,[' - Detected ' num2str(length(VariableType.aux_with_w_dependence)) ' uncertainty dependent auxilliary variables']);
    end
    if strcmp('none',ops.robust.auxreduce)
         dispsilent(ops.verbose,[' - Using possibly conservative approach to deal with uncertainty dependent auxilliary variables.']);
         dispsilent(ops.verbose,[' - (change robust.auxreduce to ''projection'', ''enumeration'' for exact solution.)'])                 
        VariableType.x_variables = unique([VariableType.aux_with_w_dependence(:)' VariableType.x_variables(:)']);
    elseif strcmp('affine',ops.robust.auxreduce)
        % Add linear feedback on all uncertainty dependent auxilliary
        % variables. The new model is dealt with as usual
        dispsilent(ops.verbose,[' - Adding affine feedback on auxilliary variables.'])
        [F_xw,xnew,info] = adjustable(F_xw,w,unique(VariableType.aux_with_w_dependence),uncertaintyGroups,ops,VariableType); 
        dispsilent(ops.verbose,[' - Feedback structure had sparsity ' num2str(info.sparsity) ' and required ' num2str(info.nvars) ' variable(s)'])
        %xnew = [];
        VariableType.x_variables = unique([VariableType.aux_with_w_dependence(:)' VariableType.x_variables(:)' getvariables(xnew)]);
    elseif strcmp('enumeration',ops.robust.auxreduce)
        % Uncertainty model is polytopic, hence enumeration can be used.
        if  isequal(Zmodel.K.s,0) & isequal(Zmodel.K.q,0)
            dispsilent(ops.verbose,[' - Using enumeration to deal with uncertainty dependent auxilliary variables']);
            ops.robust.lplp = 'enumeration';
        else
            disp([' - Cannot robustify model exactly with uncertainty dependent auxilliary variables and conic uncertainty']);
            disp([' - Change the option robust.auxreduce to ''projection'', ''affine'' or ''none''.'])
            error('robustification failed');
        end
        % User wants to do enumeration based removal of auxilliary
        % variables, hence we cannot switch later 
        ops.robust.forced_enumeration = 1;
    elseif ~(any(is(F_xw,'sdp')) | any(is(F_xw,'socp'))) | ~isempty(strfind('projection',ops.robust.auxreduce))
        if isequal(ops.robust.auxreduce,'projection')
            dispsilent(ops.verbose,[' - Projecting out auxilliary variables. This can take time...']);
        end
        F_xw = projectOut(F_xw,w,unique(VariableType.aux_with_w_dependence),uncertaintyGroups,ops);        
    else
        dispsilent(ops.verbose,[' - Cannot robustify model exactly with uncertainty dependent auxilliary variables and conic uncertainty']);
        dispsilent(ops.verbose,[' - Change the option robust.auxreduce to ''affine'' or ''none'' to compute conservative solution'])
        error('robustification failed');
    end
end

% Separate the uncertainty models accroding to uncertainty groups
separatedZmodel = separateUncertainty(Zmodel,uncertaintyGroups);

% Conversion of bounded variables that have been modelled using
% the norm operator (remove the epigraph variable to ensure explicit
% maximization is used). This will be generalized in the next version
[separatedZmodel, uncertaintyGroups] = convertUncertaintyGroups(separatedZmodel,uncertaintyGroups,VariableType);

% Code will be added to detect uncertainty cases in a more general and
% modular way. Additional code will also be added to find hidden simple
% structures, such as norm(w,1)<1, which currently is treated as a general
% polytopic uncertainty, since the expansion hides the simplicity
% 'Bounds', 'Simplex', 'Conic', 'Polytopic', '2-norm', '1-norm', 'inf-norm'
[uncertaintyTypes,separatedZmodel,uncertaintyGroups] = classifyUncertainty(separatedZmodel,uncertaintyGroups,w);

% Misplaced constraints. Really isn't uncertain, but when we expanded an
% uncertain operator, a new auxilliary variable was introduced, but has not
% made dependent on w. Hewnce, it should really be moved. Taken care
% outside now though
% if length(F_xw)>0
%     move = zeros(length(F_xw),1);
%     for i = 1:length(F_xw)
%         if all(ismember(getvariables(F_xw(i)),VariableType.x_variables))
%             move(i) = 1;
%         end
%     end
%     if any(move)
%         F_x = [F_x, F_xw(find(move))];
%         F_xw(find(move))=[];
%     end
% end
            
UncertainModel.F_x = F_x;
UncertainModel.F_xw = F_xw;
UncertainModel.h = h;
Uncertainty.F_w = F_w;
Uncertainty.Zmodel = Zmodel;
Uncertainty.separatedZmodel = separatedZmodel;
Uncertainty.uncertaintyTypes = uncertaintyTypes;
Uncertainty.separatedZmodel = separatedZmodel;
Uncertainty.uncertaintyGroups = uncertaintyGroups;

VariableType.x = recover(VariableType.x_variables);
VariableType.w = recover(VariableType.w_variables);

failure = 0;


function dispsilent(notsilent,text)
if notsilent
    disp(text);
end

function [F,w] = extractUncertain(F,w);
if isempty(w)
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = recover(getvariables(sdpvar(F(find(unc_declarations)))));
        F = F(find(~unc_declarations));
    else
        error('There is no uncertainty definition in the model.')
    end
end

function groups = findSeparate(model)

% This is an early bail-out to avoid any errors during the devlopment of
% new features. Separable constraints are only support for Polya models
if any(model.K.s > 0) %| any(model.K.q >0)
    groups{1} = 1:size(model.F_struc,2)-1;
    return
end

X = zeros(size(model.F_struc,2)-1);
top  = 1;
if model.K.f + model.K.l > 0
    A = model.F_struc(top:model.K.f+model.K.l,2:end);
    for i = 1:size(A,1)
        X(find(A(i,:)),find(A(i,:))) = 1;
    end
    top = top + model.K.f + model.K.l;
end

if any(model.K.q > 0)
    for j = 1:length(model.K.q)
        A = model.F_struc(top:top+model.K.q(j)-1,2:end);top = top + model.K.q(j);
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

if any(model.K.s > 0)
    for j = 1:length(model.K.s)
        A = model.F_struc(top:top+model.K.s(j)^2-1,2:end);top = top + model.K.s(j)^2;
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

[a,b,c,d] = dmperm(X);
for i = 1:length(d)-1
    groups{i} = sort(a(d(i):d(i+1)-1));
end

function newModels = separateUncertainty(Zmodel,uncertaintyGroups);

for i = 1:length(uncertaintyGroups)
    data = Zmodel.F_struc(:,[1 1+uncertaintyGroups{i}]);
    
    data_f = data(1:Zmodel.K.f,:);
    data_l = data(Zmodel.K.f+1:Zmodel.K.f+Zmodel.K.l,:);
    %data_q = data(Zmodel.K.f+Zmodel.K.l+1:sum(Zmodel.K.q),:);
    %data_s = data(Zmodel.K.f+Zmodel.K.l+1:sum(Zmodel.K.q)+1:end,:);
    
    eqIndex = find(any(data_f(:,2:end),2));
    liIndex = find(any(data_l(:,2:end),2));
    
    newModels{i}.F_struc = [data_f(eqIndex,:);data_l(liIndex,:)];
    newModels{i}.K.f = length(eqIndex);
    newModels{i}.K.l = length(liIndex);
    newModels{i}.K.q = [];
    top = Zmodel.K.f+Zmodel.K.l+1;
    if Zmodel.K.q(1)>0
        for j = 1:length(Zmodel.K.q)
            data_q = data(top:top+Zmodel.K.q(j)-1,:);
            if nnz(data_q(:,2:end))>0
                newModels{i}.K.q(end+1) = Zmodel.K.q(j);
                newModels{i}.F_struc = [newModels{i}.F_struc;data_q];
            end
            top = top + Zmodel.K.q(j);
        end
    end
    if isempty(newModels{i}.K.q)
        newModels{i}.K.q = 0;
    end
    newModels{i}.K.s = [];
    top = Zmodel.K.q+Zmodel.K.f+Zmodel.K.l+1;
    if Zmodel.K.s(1)>0
        for j = 1:length(Zmodel.K.s)
            data_s = data(top:top+Zmodel.K.s(j)^2-1,:);
            if nnz(data_s(:,2:end))>0
                newModels{i}.K.s(end+1) = Zmodel.K.s(j);
                newModels{i}.F_struc = [newModels{i}.F_struc;data_s];
            end
            top = top + Zmodel.K.s(j)^2;
        end
    end
    if isempty(newModels{i}.K.s)
        newModels{i}.K.s = 0;
    end
    
   % newModels{i}.K.s = 0;
    newModels{i}.variabletype = Zmodel.variabletype(uncertaintyGroups{i});
end

function Zmodel = convertuncertainty(Zmodel);
% Temoporary hack, will be generalized once the framework for multiple
% uncertainty models is supported
% We are looking for k>t, -tw<t
if size(Zmodel,1) == 1+(size(Zmodel,2)-1)*2 & Zmodel.K.f==0 & Zmodel.K.l==size(Zmodel.F_struc,1)
    n = size(Zmodel.F_struc,2)-1;
    if isequal(Zmodel.F_struc(:,2:end),sparse([zeros(1,n-1) -1;[eye(n-1);-eye(n-1)] ones(2*(n-1),1)]))
        Zmodel.F_struc = [ones(2*n,1)*Zmodel.F_struc(1,1) [eye(n);-eye(n)]];
        Zmodel.K.l = 2*n;
    end
end

function [separatedZmodel, uncertaintyGroups,VariableType] = convertUncertaintyGroups(separatedZmodel,uncertaintyGroups,VariableType);
% Temporary hack, will be generalized once the framework for multiple
% uncertainty models is supported

% We are looking for k>t, -t<w<t. This is the slightly redundant model for
for i = 1:length(separatedZmodel)
    % -k<w<k genrated when YALMIP encounters norm(w,inf)<k
    if size(separatedZmodel{i},1) == 1+(size(separatedZmodel{i},2)-1)*2 & separatedZmodel{i}.K.f==0 & separatedZmodel{i}.K.l==size(separatedZmodel{i}.F_struc,1)
        n = size(separatedZmodel{i}.F_struc,2)-1;
        if isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(1,n-1) -1;[eye(n-1);-eye(n-1)] ones(2*(n-1),1)]))
            k = separatedZmodel{i}.F_struc(1,1);
            c = separatedZmodel{i}.F_struc(2:2+n-2,1);
            separatedZmodel{i}.F_struc = [[k-c;k;c+k;k] [-eye(n);eye(n)]];
            separatedZmodel{i}.K.l = 2*n;
        end
        
    elseif separatedZmodel{i}.K.l == 1 & separatedZmodel{i}.K.q(1)>0 & length(separatedZmodel{i}.K.q(1))==1
        % Could be norm(w,2) < r.
        [n,m] = size(separatedZmodel{i}.F_struc);
        if n==m & n>=3
            n = n-2;
            
            if isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(2,n) [-1;1];eye(n) zeros(n,1)]))% isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(2,n) [-1;1];eye(n) zeros(n,1)]))
                if separatedZmodel{i}.F_struc(1,1)>0 %& nnz(separatedZmodel{i}.F_struc(2:end,1))==0
                    % the user has written norm(w-xc) < r. YALMIP will handle
                    % this using the nonlinear operator framework, and
                    % introduce a variable t and write it as norm(w)<t, t<r
                    % The variable t will be appended to the list of
                    % uncertain variables, and thus the model is a SOCP
                    % uncertainty in (w,t). We however want to interpret it
                    % as a simple quadratic model. We thus write it as
                    % (w,t)Q(w,t) < 1. Since t actually is unbounded in
                    % this form (since t is an auxilliary variable that we
                    % now through away, we arbitrarily say Q=[I/r^2 0;0 1]
                    r = separatedZmodel{i}.F_struc(1,1);
                    center = -separatedZmodel{i}.F_struc(3:end,1);
                    % separatedZmodel{i}.Q = blkdiag(r*eye(n),1);
                    % separatedZmodel{i}.center = [center;0];
                    
                    separatedZmodel{i}.K.l = 0;
                    separatedZmodel{i}.K.q = n+1;
                    % r = separatedZmodel{i}.F_struc(1,1);
                    % center = -separatedZmodel{i}.F_struc(3:end,1);
                    separatedZmodel{i}.F_struc = [[1+r^2;-2*center;1-r^2] [zeros(1,n);eye(n)*2;zeros(1,n)]];
                   % VariableType.w_variables = setdiff(VariableType.w_variables,uncertaintyGroups{1}(end));
                    uncertaintyGroups{1}(end)=[];                    
                end
            end
        end
    end
    
end


function [uncertaintyTypes,separatedZmodel,uncertaintyGroups] = classifyUncertainty(separatedZmodel,uncertaintyGroups,w)

for i = 1:length(separatedZmodel)
    uncertaintyTypes{i} = 'unclassified';
end

% Look for simplicies, which can be used in Polya
simplex_members = find_simplex_models(separatedZmodel);
for i = find(simplex_members)
    uncertaintyTypes{i} = 'simplex';
end

% Look for simple bounds, and them combine them
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if any(separatedZmodel{i}.K.q > 0) | any(separatedZmodel{i}.K.s > 0)
            simplebounds = 0;
        else
            [aux,lower,upper] = find_simple_variable_bounds(separatedZmodel{i});
            simplebounds = ~isinf(lower) & ~isinf(upper);
        end
        if all(simplebounds)
            if aux.K.l == 0
                uncertaintyTypes{i} = 'inf-norm';
                separatedZmodel{i}.lb = lower;
                separatedZmodel{i}.ub = upper;
            end
        end
    end
end

j = find(strcmp(uncertaintyTypes,'inf-norm'));
if length(j)>1
    allBounded = [];
    lb = [];
    ub = [];
    for i = 1:length(j)
        allBounded = [allBounded; uncertaintyGroups{j(i)}];
        lb = [lb;separatedZmodel{j(i)}.lb];
        ub = [ub;separatedZmodel{j(i)}.ub];
        if any(lb > ub)
            error('There are inconsistent bounds in the uncertainty model');
        end
    end
    separatedZmodel{j(1)}.lb = lb;
    separatedZmodel{j(1)}.ub = ub;
    uncertaintyGroups{j(1)} = allBounded(:)';
    keep = setdiff(1:length(separatedZmodel),j(2:end));
    separatedZmodel = {separatedZmodel{keep}};
    uncertaintyGroups = {uncertaintyGroups{keep}};
    uncertaintyTypes = {uncertaintyTypes{keep}};
end

% Look for 2-norm balls norm(w) < r, written using norm(w) or w'*w
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.s > 0) &  separatedZmodel{i}.K.l==0 & separatedZmodel{i}.K.f==0 & length(separatedZmodel{i}.K.q)==1
            % Hmm, only 1 SOCP
            if all(separatedZmodel{i}.F_struc(1,2:end)==0)
                if isequal(separatedZmodel{i}.F_struc(2:end,2:end),speye(separatedZmodel{i}.K.q(1)-1))
                    % r > norm(x)
                    uncertaintyTypes{i} = '2-norm';
                    separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end,1);
                    separatedZmodel{i}.r = sqrt(max(0,separatedZmodel{i}.F_struc(1,1)^2));
                elseif size(separatedZmodel{i}.F_struc,1)==size(separatedZmodel{i}.F_struc,2)
                    S = separatedZmodel{i}.F_struc(2:end,2:end);
                    [ii,jj,kk] = find(S);
                    if all(kk==2) & length(ii)==length(unique(ii)) & length(jj)==length(unique(jj))
                        % r > norm(x-center)
                        uncertaintyTypes{i} = '2-norm';
                        separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end,1)/2;
                        separatedZmodel{i}.r = sqrt((separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2)/4);
                    end
                elseif separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2>0
                    S = separatedZmodel{i}.F_struc(2:end,2:end);
                    [ii,jj,kk] = find(S);
                    if all(kk==2) & length(ii)==length(unique(ii)) & length(jj)==length(unique(jj))
                        uncertaintyTypes{i} = '2-norm';
                        separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end-1,1)/2;
                        separatedZmodel{i}.r = sqrt((separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2)/4);
                    end
                end
            end
        end
    end
end

% Look for quadratic constraints, other than norm models
% A bit redundant code should be integrated with the case above
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.s > 0) &  separatedZmodel{i}.K.l==0 & separatedZmodel{i}.K.f==0 & length(separatedZmodel{i}.K.q)==1
            % 1 single SOCP ||Ax+b|| <= cx+d
            A = separatedZmodel{i}.F_struc(2:end,2:end);
            b = separatedZmodel{i}.F_struc(2:end,1);
            c = separatedZmodel{i}.F_struc(1,2:end);
            d = separatedZmodel{i}.F_struc(1,1);
            if min(eig(full(A'*A-c'*c)))
                % Originates in a quadratic constraint
                rhs = c*w(uncertaintyGroups{i}(:))+d;
                lhs = A*w(uncertaintyGroups{i}(:))+b;
                uncertaintyTypes{i} = 'quadratic';                
                separatedZmodel{i}.g = rhs^2-lhs'*lhs;
                separatedZmodel{i}.center = [];
                separatedZmodel{i}.r = [];
            end
        end
    end
end

% Look for 1-norm balls |w|1 < r^2
% We are looking for the auto-generate model -t<w<t, sum(t)<s,s<r^2
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.s > 0) &  separatedZmodel{i}.K.l>0 & separatedZmodel{i}.K.f==0 & ~any(separatedZmodel{i}.K.q > 0)
            %if all(separatedZmodel{i}.F_struc(2:end,1)==0)
            if separatedZmodel{i}.F_struc(1,1)>0
                n = (size(separatedZmodel{i}.F_struc,2)-2)/2;
                if n==fix(n)
                    try
                        if all(separatedZmodel{i}.F_struc(n+2:end-1,1)==-separatedZmodel{i}.F_struc(2:2+n-1,1))
                            shouldbe = [zeros(1,n) -1 zeros(1,n);eye(n) zeros(n,1) eye(n);-eye(n) zeros(n,1) eye(n);zeros(1,n) 1 -ones(1,n)];%; zeros(n,n+1)  eye(n)];;zeros(1,n) 1 zeros(1,n)];
                            if  isequal(full(separatedZmodel{i}.F_struc(:,2:end)),shouldbe)
                                %if all(separatedZmodel{i}.F_struc(2:end,1)==0)
                                uncertaintyTypes{i} = '1-norm';
                                separatedZmodel{i}.r = separatedZmodel{i}.F_struc(1,1);
                                separatedZmodel{i}.center =separatedZmodel{i}.F_struc(n+2:end-1,1);
                                %end
                            end
                        end
                    catch
                        %FIX ME, caught by -1<w<1, sum(w)<1
                    end
                end
            end
            %end
        end
    end
end


% Look for polytopes
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.q > 0) | any(separatedZmodel{i}.K.s > 0)
        end
    end
end

function Fnew = projectOut(F,w,newAuxVariables,uncertaintyGroups,ops)

w_variables = getvariables(w);

Graph = yalmip('getdependence');
for i = 1:length(uncertaintyGroups)
    Graph(w_variables(uncertaintyGroups{i}),w_variables(uncertaintyGroups{i}))=1;
end
aux_and_w = union(newAuxVariables,w_variables);
F_lp = F(find(is(F,'elementwise')));
X = sdpvar(F_lp);
Xvars = getvariables(X);
G = getbase(X);G = G(:,2:end);
for i = 1:size(G,1)
    j = Xvars(find(G(i,:)));
    Graph(j,j) = 1;
end
Graph(:,setdiff(1:size(Graph,2),aux_and_w))=0;
Graph(setdiff(1:size(Graph,2),aux_and_w),:)=0;
Graph2 = Graph(1:max(aux_and_w),1:max(aux_and_w));
Graph2 = Graph2 + speye(length(Graph2));
[aa,bb,cc,dd] = dmperm(Graph2);
commonProjections = {};
for r = 1:length(dd)-1
    comps = dd(r):dd(r+1)-1;
    comps = intersect(aa(comps),newAuxVariables);
    if ~isempty(comps)
        commonProjections{end+1}  = comps;
    end
end

if ops.verbose
    disp([' - * Partitioned these into ' num2str(length(commonProjections)) ' group(s)']);
end

keep = ones(length(F),1);
Fnew = [];
started = 0;
for clique = 1:length(commonProjections)
    F_lp = [];
    for i = 1:length(F)
        F_vars = getvariables(F(i));
        var_w = intersect(F_vars,w_variables);
        var_aux = intersect(F_vars,commonProjections{clique});
        
        if is(F(i),'elementwise')
            if any(ismember(getvariables(F(i)),commonProjections{clique}))
                F_lp = [F_lp, F(i)];
                keep(i) = 0;
                if ~started
                    started = 1;
                    if ops.verbose
                        disp([' - * Performing projections of uncertain graph variables...']);
                    end
                end
            end
        end
    end
    if ~isempty(F_lp)
        X = sdpvar(F_lp);
        [Ab] = getbase(X);
        b = Ab(:,1);
        A = -Ab(:,2:end);
        allVariables = getvariables(X);
        newAuxVariables = intersect(commonProjections{clique},allVariables);
        j = [];
        for i = 1:length(newAuxVariables)
            j(i) = find(allVariables==newAuxVariables(i));
        end
        
        [A,b] = fourierMotzkin(A,b,sort(j));
        A(:,j)=[];
        left = recover(setdiff(allVariables,newAuxVariables));
        X = b-A*left(:);
        Fnew = [Fnew,[[X >= 0] : 'Projected uncertain']];
    end
end
if started & ops.verbose
    d = length(sdpvar(Fnew))-length(sdpvar(F_lp));
    if d>0
    disp([' - * Done with projections (generated ' num2str(d) ' new constraints)']);
    elseif d<0
        disp([' - * Done with projections (actually reduced model with ' num2str(-d) ' constraints)']);
    else
        disp([' - * Done with projections (model size unchanged)']);
    end
end
Fnew = [Fnew,F(find(keep))];

function [A,b] = fourierMotzkin(A,b,dim)
% Brute-force projection through Fourier-Motzkin
for i = dim(:)'
    
    a = A(:,i);
    zero = find(a == 0);
    pos = find(a > 0);
    neg = find(a < 0);
    
    T = spalloc(length(zero) + length(pos)*length(neg),size(A,1),length(zero)+length(pos)*length(neg)*2);
    
    row = 1;
    for j = zero(:)'
        T(row,j) = 1;
        row = row+1;
    end
    for j = pos(:)'
        for k = neg(:)'
            T(row,j) = -a(k);
            T(row,k) = a(j);
            row = row+1;
        end
    end    
    A = T*A;
    b = T*b;
end

function [Fnew,xnew,info] = adjustable(F,w,newAuxVariables,uncertaintyGroups,ops,VariableType)

L = sdpvar(length(recover(newAuxVariables)),length(w),'full');
L = L.*VariableType.Graph(VariableType.aux_with_w_dependence,VariableType.w_variables);
y0 = sdpvar(length(recover(newAuxVariables)),1);
Fnew = replace(F,recover(newAuxVariables),y0+L*w);
xnew = recover([getvariables(y0) getvariables(L)]);
info.nvars = length(getvariables(L));
info.sparsity = length(getvariables(L))/prod(size(L));

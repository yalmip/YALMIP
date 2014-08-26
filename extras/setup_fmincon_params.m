function model = setup_fmincon_params(model)

monomtable = model.monomtable;
nonlinearindicies = model.nonlinearindicies;
linearindicies = model.linearindicies;

model.evalinobjective  = 0;

if ~isempty(model.evalMap)
    [i,j,k] = find(model.monomtable(:,model.evalVariables));
    evalInvolved = union(model.evalVariables,i);
    model.evalinobjective = any(model.c(evalInvolved)) | any(model.Q(evalInvolved,1));
    if size(model.Anonlineq,1)>0
        if nnz(model.Anonlineq(:,evalInvolved))>0
            model.evalinconstraint = 1;
        end
    end
    if size(model.Anonlinineq,1)>0
        if nnz(model.Anonlinineq(:,evalInvolved))>0
            model.evalinconstraint = 1;
        end
    end
    if any(model.K.q)
        if nnz(model.F_struc(1+model.K.f + model.K.l:end,1+evalInvolved)) > 0
            model.evalinconstraint = 1;
        end
    end
    % Speed up callbacks by removing last marker argument (only used in
    % expandmodel etc) and figuring out what R^m -> R^n operators computes
    for i = 1:length(model.evalMap)        
        model.evalMap{i}.prearg = {model.evalMap{i}.fcn,model.evalMap{i}.arg{1:end-1}};
        model.evalMap{i}.prearg{1+model.evalMap{i}.argumentIndex} = [];
        if ~isfield(model.evalMap{i},'computes')
            model.evalMap{i}.computes = model.evalVariables(i);
        end
    end
end

% Figure out if YALMIP easily can compute the gradient of the objective
% This will done completely general later
model.SimpleLinearObjective = 0;    % Linear
model.SimpleQuadraticObjective = 0; % Quadratic
model.SimpleNonlinearObjective = 1; % Polynomial
if isempty(model.evalMap)
    if nnz(model.c(nonlinearindicies)) == 0
        if (nnz(model.Q)==0)
            model.SimpleLinearObjective = 1;
        else
            if nnz(model.Q(nonlinearindicies,nonlinearindicies))==0
                model.SimpleQuadraticObjective = 1;
            end
        end
    end
elseif ~model.evalinobjective
    if nnz(model.c(nonlinearindicies)) == 0
        if (nnz(model.Q)==0)
            model.SimpleLinearObjective = 1;
        else
            if nnz(model.Q(nonlinearindicies,nonlinearindicies))==0
                model.SimpleQuadraticObjective = 1;
            end
        end
    elseif nnz(model.Q)==0
        r = find(model.c);
        if all(model.variabletype(r)<=2)
            if isempty(intersect(r,model.evalVariables))
                % Simple quadratic function
                for i = r(:)'
                    if model.variabletype(i)==2
                        j = find(model.monomtable(i,:));
                        model.Q(j,j) = model.c(i);
                        model.c(i) = 0;
                    elseif model.variabletype(i)==1
                        j = find(model.monomtable(i,:));
                        model.Q(j(1),j(2)) = model.c(i)/2;
                        model.Q(j(2),j(1)) = model.c(i)/2;
                        model.c(i) = 0;
                    end
                end
                model.SimpleQuadraticObjective = 1;
            end
        end
    end
else
    model.SimpleNonlinearObjective = 0;
end

model.linearconstraints = isempty(model.Anonlinineq) & isempty(model.Anonlineq) & nnz(model.K.q)==0;
model.nonlinearinequalities = ~isempty(model.Anonlinineq);
model.nonlinearequalities = ~isempty(model.Anonlineq);
if any(model.K.q)
    if nnz(model.F_struc(1+model.K.f + model.K.l:end,1+model.nonlinearindicies)) > 0
        model.nonlinearcones = 1;
    else
        model.nonlinearcones = 0;
    end
else
    model.nonlinearcones = 0;
end

% Structure for fast evaluation of monomial terms in differentiation
 if isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities | model.nonlinearcones) & ~isfield(model,'fastdiff')     
    allA = [model.Anonlineq;model.Anonlinineq];
    dgAll = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    
    allDerivemt = [];
    news = [];
    c = [];
    mtNonlinear = mtNonlinear(:,linearindicies);
    mtNonlinear = mtNonlinear';
    
    [jj,ii,val] = find(mtNonlinear');
    for k = 1:length(ii)
        i = ii(k);
        j = jj(k);
        s=mtNonlinear(:,j);
        c = [c;s((i))];
        s((i)) = s((i))-1;
        allDerivemt = [allDerivemt s(:)];
        news = [news;j i];
    end
    
    allDerivemt = allDerivemt';
    
    model.fastdiff.news = news;
    model.fastdiff.allDerivemt = allDerivemt;
    model.fastdiff.c = c;
    model.fastdiff.univariateDifferentiates = 0;
        
    a1 =  model.fastdiff.news(1:length(model.fastdiff.c),2); 
    a2 =  model.nonlinearindicies(model.fastdiff.news(1:length(model.fastdiff.c),1))'; 
    zzz = [ones(length(a1),1)]; 
    nn = max(max(length(model.linearindicies)),max(news(:,2))); 
    mm = max(max(linearindicies),max(model.nonlinearindicies(news(:,1))));
    a1f = [a1(:);(1:length(model.linearindicies))']; 
    a2f = [a2(:);model.linearindicies(:)];    
    zzzf = [zzz;ones(length(linearindicies),1)]; 
    
    %model.fastdiff.newdxx = sparse(a1f,a2f,zzzf,nn,mm); 
    %model.fastdiff.linear_in_newdxx = sub2ind([nn mm],a1(:),a2(:));
    model.fastdiff.newdxx = sparse(a1f,a2f,zzzf,nn,mm)'; 
    model.fastdiff.linear_in_newdxx = sub2ind([mm nn],a2(:),a1(:));
    model.fastdiff.newdxx(model.fastdiff.linear_in_newdxx)=1;
    
    if all(sum(allDerivemt | allDerivemt,2)==1)
        model.fastdiff.univariateDifferentiates = 1;
        [i,j,k] = find(allDerivemt');
        if any(diff(j)<0)
            error;
        end
        model.fastdiff.univariateDiffMonom = i(:);
        model.fastdiff.univariateDiffPower = k(:);
    end
 else
      allA = [model.Anonlineq;model.Anonlinineq]; 
      if any(model.K.q)
          allA = [allA;model.F_struc(1+model.K.f + model.K.f:end,2:end)];
      end
      requested = any(allA',2); 
      [i,j,k] = find((model.deppattern(find(requested),:))); 
      requested(j) = 1; 
      model.fastdiff.requested = requested;     
 end
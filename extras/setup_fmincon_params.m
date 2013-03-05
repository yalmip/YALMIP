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
    % Speed up callbacks by removing last marker argument (only used in
    % expandmodel etc) and figuring out what R^m -> R^n operators computes
    for i = 1:length(model.evalMap)
        %model.evalMap{i}.prearg = {model.evalMap{i}.fcn,[],model.evalMap{i}.arg{2:end-1}};
        model.evalMap{i}.prearg = {model.evalMap{i}.fcn,model.evalMap{i}.arg{1:end-1}};
        model.evalMap{i}.prearg{1+model.evalMap{i}.argumentIndex} = [];
        if isfield(model.evalMap{i},'computes')
%             temp = zeros(1,length(model.evalMap{i}.computes));
%             for j = 1:length(model.evalMap{i}.computes)
%                 temp(j) = find(model.evalMap{i}.computes(j) == model.used_variables);
%             end
%             model.evalMap{i}.computes = temp;
        else
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

model.linearconstraints = isempty(model.Anonlinineq) & isempty(model.Anonlineq);
model.nonlinearinequalities = ~isempty(model.Anonlinineq);
model.nonlinearequalities = ~isempty(model.Anonlineq);

% Structure for fast evaluation of monomial terms in differentiation
 if isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities) & ~isfield(model,'fastdiff')     
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

 
%     
%     for i = 1:length(linearindicies)
%         for j = 1:size(mtNonlinear,2)
%             if mtNonlinear(i,j);
%                 s=mtNonlinear(:,j);               
%                 c = [c;s((i))];               
%                 s((i)) = s((i))-1;                
%                 allDerivemt = [allDerivemt s(:)];
%                 news = [news;j i];
%             end
%         end
%     end
    
    
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
    model.fastdiff.newdxx = sparse(a1f,a2f,zzzf,nn,mm); 
    model.fastdiff.linear_in_newdxx = sub2ind([nn mm],a1(:),a2(:));
    
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
      requested = any(allA',2); 
      [i,j,k] = find((model.deppattern(find(requested),:))); 
      requested(j) = 1; 
      model.fastdiff.requested = requested;     
 end
    
    





































% % 
% % 
% % 
% % function constant_data = setup_fmincon_params(model)
% % 
% % monomtable = model.monomtable;
% % nonlinearindicies = model.nonlinearindicies;
% % linearindicies = model.linearindicies;
% % 
% % constant_data = model;
% % 
% % constant_data.nonlinearineval = 0;
% % % Check if there are any nonlinear expression in evaluation based operators
% % if ~isempty(model.evalMap)
% %     temp = [];
% %     for i = 1:length(model.evalMap)
% %         temp = [temp model.evalMap{i}.variableIndex];
% %     end
% %     if any(model.variabletype(temp))
% %         constant_data.nonlinearineval = 1;
% %     end
% % end
% % constant_data.evalinconstraint = 0;
% % constant_data.evalinobjective = 0;
% % if ~isempty(model.evalMap)
% %     [i,j,k] = find(model.monomtable(:,model.evalVariables));
% %     evalInvolved = union(model.evalVariables,i);
% %     if size(model.Anonlinineq,1)>0
% %         if nnz(model.Anonlinineq(:,evalInvolved))>0
% %             constant_data.evalinconstraint = 1;
% %         end
% %     end
% %     constant_data.evalinobjective = any(model.c(evalInvolved)) | any(model.Q(evalInvolved,1));
% %     if size(model.Anonlineq,1)>0
% %         if nnz(model.Anonlineq(:,evalInvolved))>0
% %             constant_data.evalinconstraint = 1;
% %         end
% %     end
% %     for i = 1:length(model.evalMap)
% %         constant_data.evalMap{i}.prearg = {model.evalMap{i}.fcn,[],model.evalMap{i}.arg{2:end-1}};
% %     end
% % end
% % constant_data.monominobjective = 0;
% % if any(model.variabletype)
% %     monoms = find(model.variabletype);
% %     if any(model.c(monoms))
% %         % Monomials in objective
% %         constant_data.monominobjective = 1;
% %         %    else
% %         % Maybe indirectly via evaluation variables
% %         %        [i,j,k] = find(model.monomtable(:,model.evalVariables));
% %     end
% % end
% % constant_data.monominconstraint = 0;
% % if any(model.variabletype)
% %     monoms = find(model.variabletype);
% %     if ~isempty(monoms)
% %         AA = [model.Anonlinineq;model.Anonlineq];
% %         if ~isempty(AA)
% %             if nnz(AA(:,monoms))>0
% %                 % Monomials in objective
% %                 constant_data.monominconstraint = 1;
% %                 %    else
% %                 % Maybe indirectly via evaluation variables
% %                 %        [i,j,k] = find(model.monomtable(:,model.evalVariables));
% %             end
% %         end
% %     end
% % end
% % 
% % % Figure out if YALMIP easily can compute the gradient of the objective
% % % This will done completely general later
% % constant_data.SimpleLinearObjective = 0;
% % constant_data.SimpleQuadraticObjective = 0;
% % constant_data.SimpleNonlinearObjective = 1;
% % constant_data.SimpleNonlinearConstraints = 0;
% % if isempty(model.evalMap)
% %     if nnz(model.c(nonlinearindicies)) == 0
% %         if (nnz(model.Q)==0)
% %             constant_data.SimpleLinearObjective = 1;
% %         else
% %             if nnz(model.Q(nonlinearindicies,nonlinearindicies))==0
% %                 constant_data.SimpleQuadraticObjective = 1;
% %             end
% %         end
% %     end
% %     if isequal(model.K.s,0) & isequal(model.K.q,0) & isequal(model.K.r,0)
% %         constant_data.SimpleNonlinearConstraints = 1;
% %     end
% % elseif ~constant_data.evalinobjective
% %     if nnz(model.c(nonlinearindicies)) == 0
% %         if (nnz(model.Q)==0)
% %             constant_data.SimpleLinearObjective = 1;
% %         else
% %             if nnz(model.Q(nonlinearindicies,nonlinearindicies))==0
% %                 constant_data.SimpleQuadraticObjective = 1;
% %             end
% %         end
% %     end
% % else
% %     constant_data.SimpleNonlinearObjective = 0;
% % end
% % 
% % if all(model.variabletype <= 2)
% %     aux1 = compile_nonlinear_table(model);
% %     constant_data.bilinears = aux1.bilinears;
% %     %    [aux1,constant_data.bilinears,aux2] = compile_nonlinear_table(model);
% % else
% %     constant_data.bilinears = [];
% % end
% % 
% % constant_data.linearconstraints = isempty(constant_data.Anonlinineq) & isempty(constant_data.Anonlineq) & isequal(constant_data.K.q,0) & isequal(constant_data.K.s,0);
% % constant_data.nonlinearinequalities = ~isempty(constant_data.Anonlinineq);
% % constant_data.nonlinearequalities = ~isempty(constant_data.Anonlineq);


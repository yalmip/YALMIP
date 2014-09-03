function [VariableType,F_x,F_w,F_xw,h] = robust_classify_variables_newest(F,h,ops,w);

Dependency = iterateDependance( yalmip('monomtable') | yalmip('getdependence') |  yalmip('getdependenceUser'));
DependsOnw = find(any((Dependency(:,getvariables(w))),2));

h_variables = getvariables(h);
h_w = find(ismember(h_variables,DependsOnw));
if ~isempty(h_w)
    base = getbase(h);
    h0 = base(1);
    base = base(2:end);base = base(:);
    sdpvar t
    
    F = [F,base(h_w(:))'*recover(h_variables(h_w)) <= t];
    base(h_w) = 0;
    h = base(:)'*recover(h_variables) + t;
    
    Dependency = iterateDependance(yalmip('monomtable') | yalmip('getdependence') | yalmip('getdependenceUser'));
    DependsOnw = find(any((Dependency(:,getvariables(w))),2));    
end

DoesNotDependOnw = find(~any((Dependency(:,getvariables(w))),2));
[notused,x_variables] = find(Dependency(DoesNotDependOnw,:));

F_w = [];
F_x = [];
F_xw = [];
for i = 1:length(F)
    F_vars = getvariables(F(i));
    F_vars = find(any((Dependency(F_vars,:)),1));
    if all(ismember(F_vars,DependsOnw))
        F_w = F_w + F(i);
    elseif all(ismember(F_vars,DoesNotDependOnw))
       F_x = F_x + F(i);
    else
       F_xw = F_xw + F(i);
    end            
end
ops.removeequalities = 0;
[F_x,failure,cause] = expandmodel(F_x,h,ops);
[F_w,failure,cause] = expandmodel(F_w,[],ops);
ops.expandbilinear = 1;
ops.reusemodel = 1; % Might be case x+norm(w)<1, norm(w)<1
[F_xw,failure,cause] = expandmodel(F_xw,h,ops,w);

w_variables = depends(F_w);
x_variables = unique([depends(F_x) depends(F_xw) depends(h)]);
x_variables = setdiff(x_variables,w_variables);

% After exanding the conic represntable, we have introduced new variables
Dependency = iterateDependance(yalmip('monomtable') | yalmip('getdependence') | yalmip('getdependenceUser'));

auxiliary = unique([yalmip('extvariables') yalmip('auxvariables')]);
if 1%~isempty(auxiliary)
    DependsOnw = find(any((Dependency(:,getvariables(w))),2));
    DoesNotDependOnw = find(~any((Dependency(:,getvariables(w))),2));
    
    temp = intersect(DependsOnw,x_variables);
    x_variables = setdiff(x_variables,DependsOnw);
    aux_with_w_dependence = temp;  
else
    aux_with_w_dependence = [];   
end

% aux_w_or_w = union(aux_with_w_dependence,w_variables);
% old_w_variables = [];
% while ~isequal(w_variables,old_w_variables);
%     old_w_variables = w_variables;
%     for i = 1:length(F_xw)
%         if all(ismember(depends(F_xw(i)),aux_w_or_w))
%             if ~any(ismember(depends(F_xw(i)),x_variables))
%                 new_w = intersect(depends(F_xw(i)),aux_w_or_w);
%                 w_variables = union(w_variables,new_w);
%                 aux_with_w_dependence = setdiff(aux_with_w_dependence,new_w);
%                 goon = 1;
%             end
%         end
%     end
%     aux_w_or_w = union(aux_with_w_dependence,w_variables);  
% end    

x = recover(x_variables);
w = recover(w_variables);

if ~isempty(F_xw)
    F_xw_scalar =  F_xw(find(is(F_xw,'elementwise') | is(F_xw,'equality')));
    F_xw_multi = F_xw - F_xw_scalar;
else
    F_xw_scalar = [];
    F_xw_multi = F_xw;
end

[MonomTable,Nonlinear] = yalmip('monomtable');
Dependency = yalmip('getdependenceUser');
evar = yalmip('extvariables');
if length(F_xw_scalar)>0
    % Optimize dependency graph
    X = sdpvar(F_xw_scalar);
    Xvar = getvariables(X);
    Xbase = getbase(X);Xbase = Xbase(:,2:end);
    for i = 1:size(Xbase,1)
        used = Xvar(find(Xbase(i,:)));
        if any(Nonlinear(used))
        used = find(any(MonomTable(used,:),1));
        end
        auxUsed = intersect(used,aux_with_w_dependence);
      
        if ~isempty(auxUsed)
            wUsed = intersect(used,w_variables);
            if ~isempty(wUsed)
                Dependency(auxUsed,wUsed) = 1;
            end       
            eUsed = intersect(used,evar);
            if ~isempty(eUsed)
                Dependency(eUsed,auxUsed) = 1;
            end
        end
    end      
end
if length(F_xw_multi) > 0
    for i = 1:length(F_xw_multi)
        used = getvariables(F_xw_multi(i));
        used = find(any(MonomTable(used,:),1));
        auxUsed = intersect((used),aux_with_w_dependence);
        wUsed = intersect((used),w_variables);
        if ~isempty(auxUsed) & ~isempty(wUsed)
            Dependency(auxUsed,wUsed) = 1;
        end
        eUsed = intersect((used),evar);
        if ~isempty(auxUsed) & ~isempty(eUsed)
            Dependency(eUsed,auxUsed) = 1;
        end
    end
end
UserDependency = yalmip('getdependenceUser');
fixed = find(any(UserDependency,2));
Dependency(fixed,:) = UserDependency(fixed,:);
Dependency = iterateDependance(Dependency);

VariableType.Graph = Dependency;
VariableType.x_variables = x_variables;
VariableType.w_variables = w_variables;
VariableType.aux_with_w_dependence = aux_with_w_dependence;

function Graph = iterateDependance(Graph)

Graph = Graph + speye(length(Graph));
Graph0 = double(Graph*Graph ~=0);
while ~isequal(Graph,Graph0)
    Graph = Graph0;
    Graph0 = double(Graph*Graph~=0);  
end

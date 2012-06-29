function simplex_model = find_simplex_models(p);

for i = 1:length(p)
    simplex_model(i)= 0;
    if p{i}.K.f == 0
        continue
    elseif any(p{i}.K.q > 0) | any(p{i}.K.s > 0)
        continue
    else
        aux = p{i};
        aux.F_struc(1:p{i}.K.f,:) = [];
        aux.K.f = 0;
        [aux,lower,upper] = find_simple_variable_bounds(aux);
        if all(lower == 0) | all(upper == 1) & p{i}.K.f == 1 & aux.K.l == 0
            simplex_model(i)=1;
        end
    end
end

% 
% 
% function [p,simplex_model] = find_simplex_models(p);

% simplex_model = [];
% if p.K.f == 0
%     return
% elseif any(p.K.q > 0) | any(p.K.s > 0)
%     return
% else
%     aux = p;
%     aux.F_struc(1:p.K.f,:) = [];
%     aux.K.f = 0;
%     [aux,lower,upper] = find_simple_variable_bounds(aux);
%     if all(lower == 0) | all(upper == 1)
%         % b = A*x
%         b = p.F_struc(1:p.K.f,1);
%         A = -p.F_struc(1:p.K.f,2:end);
%         for i = 1:length(b)
%             if abs(b)==1
%                 if all(A(i,:) == 1/abs(b(i)))
%                     simplex_model{end+1} = find(A(i,:));
%                 end
%             end
%         end
%     end
%     return
% end


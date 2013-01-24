function [NegativeSemiVar,H,c,A,lb,ub,semicont_variables] = negateNegativeSemiContVariables(H,c,A,lb,ub,semicont_variables,Qi);
%function [NegativeSemiVar,H,c,Aineq,Aeq,lb,ub,semicont_variables] = negateNegativeSemiContVariables(H,c,Aineq,Aeq,lb,ub,semicont_variables,Qi);

NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(lb(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = ub(semicont_variables(NegativeSemiVar));
        ub(semicont_variables(NegativeSemiVar)) = -lb(semicont_variables(NegativeSemiVar));
        lb(semicont_variables(NegativeSemiVar)) = -temp;
        if ~isempty(A)
             A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        end
%         if ~isempty(Aineq)
%             Aineq(:,semicont_variables(NegativeSemiVar)) = -Aineq(:,semicont_variables(NegativeSemiVar));
%         end
%         if ~isempty(Aeq)
%             Aeq(:,semicont_variables(NegativeSemiVar)) = -Aeq(:,semicont_variables(NegativeSemiVar));
%         end
%         c(semicont_variables(NegativeSemiVar)) = -c(semicont_variables(NegativeSemiVar));
        if ~isempty(H)
            H(:,semicont_variables(NegativeSemiVar)) = -H(:,semicont_variables(NegativeSemiVar));
            H(semicont_variables(NegativeSemiVar),:) = -H(semicont_variables(NegativeSemiVar),:);
        end
        
        for i = 1:length(Qi)
            Qi{i}(:,semicont_variables(NegativeSemiVar)) = -Qi{i}(:,semicont_variables(NegativeSemiVar));
            Qi{i}(semicont_variables(NegativeSemiVar),:) = -Qi{i}(semicont_variables(NegativeSemiVar),:);
        end
    end
end

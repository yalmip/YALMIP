function x0 = propagatex0(x0,Ab,K);
if size(Ab,1) > 0 & ~all(isnan(x0)) & any(isnan(x0))    
    % Only some of the initial values have been set. Try to complete them
    % by looking at equalities  Ax == b
    A = Ab(1:K.f,2:end);
    b = -Ab(1:K.f,1);
    x0 = x0(:);   
    for i = 1:K.f
        a = A(i,:);
        [~,j,s] = find(a);        
        notSet = find(isnan(x0(j)));
        if length(notSet) == 1
            j = j(:);
            theRest = find(~isnan(x0(j)));
           % a(notset)*x0(notset) + a(rest)*x0(rest) == b
            x0(j(notSet)) = (b(i) - a(j(theRest))*x0(j(theRest)))/a(j(notSet));
        end
    end   
end
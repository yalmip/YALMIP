function Z = nonlinearreplace(X,Y,W)

% Very slow, but easily coded. Hopefully not used much

% Special case for replace(x,y,-y)
if isequal(getbase(Y),[0 1]) && isequal(getbase(W),[0 -1]) && isequal(getvariables(Y),getvariables(W)) && is(Y,'linear')
    mt = yalmip('monomtable');    
    xv = getvariables(X);
    yv = getvariables(Y);
    if all(mt(:,yv) == fix(mt(:,yv)) )
        signs = ones(length(xv),1);
        for i = 1:length(xv)
            pow = mt(xv(i),yv);
            if ~even(pow)
                signs(i) = -1;
            end        
        end
        Z = X;
        Z = setbase(X,getbase(X)*diag([1;signs]));
        return
    end
end

U = sdpvar(length(W),1);
for kk = 1:length(Y)
    Z = [];
    for ii = 1:size(X,1)
        temp = [];
        for jj = 1:size(X,2);

            [coeffs,base] = coefficients(X(ii,jj),Y(kk));

            newp = 0;
            for i = 1:length(base)
                newp = newp + coeffs(i)*U(kk)^degree(base(i));
            end

            temp = [temp newp];
        end
        Z = [Z;temp];
    end
    X = Z;
end

Y = U;
for kk = 1:length(Y)
    Z = [];
    for ii = 1:size(X,1)
        temp = [];
        for jj = 1:size(X,2);

            [coeffs,base] = coefficients(X(ii,jj),Y(kk));

            newp = 0;
            for i = 1:length(base)
                newp = newp + coeffs(i)*W(kk)^degree(base(i));
            end

            temp = [temp newp];
        end
        Z = [Z;temp];
    end
    X = Z;
end




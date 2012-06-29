function Z = nonlinearreplace(X,Y,W)

% Very slow, but easily coded. Hopefully not used much
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




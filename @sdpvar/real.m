function X = real(X)
%REAL (overloaded)

X.basis = real(X.basis);
X = clean(X);
if isa(X,'sdpvar')
    X.conicinfo = [0 0];
    if length(X.midfactors)>0
        newleftfactors = {};
        newmidfactors = {};
        newrightfactors = {};
        for i = 1:length(X.midfactors)
            ar = real(X.leftfactors{i});
            ai = imag(X.leftfactors{i});
            br = real(X.rightfactors{i});
            bi = imag(X.rightfactors{i});
            
            if nnz(ar)>0 & nnz(br)>0
                newleftfactors{end+1} = ar;
                newmidfactors{end+1} = X.midfactors{i};
                newrightfactors{end+1} = br;
            end
            
            if nnz(ai)>0 & nnz(bi)>0
                newleftfactors{end+1} = -ai;
                newmidfactors{end+1} = X.midfactors{i};
                newrightfactors{end+1} = bi;
            end
            
        end
        X.leftfactors = newleftfactors;
        X.midfactors = newmidfactors;
        X.rightfactors = newrightfactors;
    end
end
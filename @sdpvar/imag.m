function X = imag(X)
%IMAG (overloaded)

X.basis = imag(X.basis);
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
            
            if nnz(ai)>0 & nnz(br)>0
                newleftfactors{end+1} = ai;
                newmidfactors{end+1} = X.midfactors{i};
                newrightfactors{end+1} = br;
            end
            
            if nnz(ar)>0 & nnz(bi)>0
                newleftfactors{end+1} = ar;
                newmidfactors{end+1} = X.midfactors{i};
                newrightfactors{end+1} = bi;
            end
            
        end
        X.leftfactors = newleftfactors;
        X.midfactors = newmidfactors;
        X.rightfactors = newrightfactors;
    end
end

%    
% Y = X;
% x_lmi_variables = X.lmi_variables;
% lmi_variables = [];
% n = X.n;
% m = X.m;
% imagX = imag(X.basis(:,1));
% Y.basis = imagX(:);
% 
% j = 1;
% for i = 1:length(x_lmi_variables)
%     imagX = imag(X.basis(:,i+1));
%     if (norm(imagX,inf)>0)
%         Y.basis(:,j+1) = imagX(:);
%         lmi_variables = [lmi_variables x_lmi_variables(i)];
%         j = j+1;
%     end
% end
% if isempty(lmi_variables)
%     Y = full(reshape(Y.basis,n,m));
% else
%     Y.lmi_variables = lmi_variables;
%     % Reset info about conic terms
%     Y.conicinfo = [0 0];
% end
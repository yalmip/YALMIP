function X = imag(X)
%IMAG (overloaded)

% Author Johan Löfberg
% $Id: imag.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

X.basis = imag(X.basis);
X = clean(X);
if isa(X,'sdpvar')
   X.conicinfo = [0 0];
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
function X = real(X)
%REAL (overloaded)

% Author Johan Löfberg 
% $Id: real.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   

%Y = X;
%x_lmi_variables = X.lmi_variables;
%lmi_variables = [];
%n = X.n;
%m = X.m;

X.basis = real(X.basis);
X = clean(X);
if isa(X,'sdpvar')
   X.conicinfo = [0 0];
end
   
% realX = real(X.basis(:,1));
% Y.basis = realX(:);
% j = 1;
% for i = 1:length(x_lmi_variables)
% 	realX = real(X.basis(:,i+1));
% 	if (norm(realX,inf)>0)
% 		Y.basis(:,j+1) = realX(:);
% 		lmi_variables = [lmi_variables x_lmi_variables(i)];
% 		j = j+1;
% 	end
% end   
% 
% if isempty(lmi_variables)
%     Y = full(reshape(Y.basis,n,m));
% else
%     Y.lmi_variables = lmi_variables;
%     % Reset info about conic terms
%     Y.conicinfo = [0 0];
% end

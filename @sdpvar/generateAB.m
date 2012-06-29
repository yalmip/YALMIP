function [A,b,lmi_variables] = generateAB(varargin)
%GENERATEAB Internal function to generates linear equation system

% Author Johan Löfberg
% $Id: generateAB.m,v 1.4 2006-07-26 20:17:58 joloef Exp $

if mod(nargin,2)~=0
    error('Internal error in generateAB. Please report.');
end

lmi_variables = [];
n_eq = 0;
for r = 1:nargin/2
    [n,m] = size(varargin{1+2*(r-1)});
    n_eq = n_eq + n*m;
    lmi_variables = [lmi_variables varargin{1+2*(r-1)}.lmi_variables(:)'];
end
lmi_variables = uniquestripped(lmi_variables);

A = spalloc(n_eq,length(lmi_variables),n_eq);
b = spalloc(n_eq,1,n_eq);

btop = 1;
atop = 1;
for r = 1:nargin/2
    X  = varargin{2*r-1};
    [n,m] = size(X);
    used_variables = find(ismembc(lmi_variables,X.lmi_variables));
    A(atop:atop+n*m-1,used_variables) = X.basis(:,2:end);
    b(btop:btop+n*m-1) =varargin{2*r}(:)-X.basis(:,1);
    atop = atop+n*m;
    btop = btop+n*m;
end

function Y=sdpvarfun(varargin)
%SDPVARFUN Applies operator on matrix variable
%
% Y = sdpvarfun('op',X)  Applies the Matlab operator 'op' on X
%
% Example : Y = sdpvarfun('triu',X)  Extracts upper triangluar part

op = varargin{1};
if ~isa(op,'char')
  error('First argument must be a string with the function name')
end

X = varargin{2};
Y = X;

x_lmi_variables = X.lmi_variables;
lmi_variables = [];

Y.basis = [];
opX = feval(op,reshape(X.basis(:,1),X.dim(1),X.dim(2)));
Y.basis = opX(:);

j = 1;
for i = 1:length(x_lmi_variables)
  opX = feval(op,reshape(X.basis(:,i+1),X.dim(1),X.dim(2)));
  if (norm(opX,inf)>0)
    Y.basis(:,j+1) = opX(:);
    lmi_variables = [lmi_variables x_lmi_variables(i)];
    j = j+1;
  end
end  
Y.dim(1) = size(opX,1);
Y.dim(2) = size(opX,2);
Y.lmi_variables = lmi_variables;
% Reset info about conic terms
Y.conicinfo = [0 0];

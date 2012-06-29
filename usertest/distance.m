function dist=distance(X,x);
% dist=distance(X,x)
%
% computes the pairwise squared distance matrix between any column vectors in X and
% in x
%
% INPUT:
%
% X     dxN matrix consisting of N column vectors
% x     dxn matrix consisting of n column vectors
%
% OUTPUT:
%
% dist  Nxn matrix 
%
% Example:
% Dist=distance(X,X);
% is equivalent to
% Dist=distance(X);
%
[D,N] = size(X);
     
try 
 if(nargin>=2)
    % PAIRWISE DISTANCES
  [D,n] = size(x);
  X2 = sum(X.^2,1);
  x2 = sum(x.^2,1);
  dotProd = X'*x; 
  dist = repmat(x2,N,1)+repmat(X2',1,n)-2*dotProd;
 else    
  [D,N] = size(X);
  if(exist('addv') & exist('addh'))
   X2 = sum(X.^2,1);
   dist=addh(addv(-2*X'*X,X2),X2);
  else
   X2 = repmat(sum(X.^2,1),N,1);
   dist = X2+X2'-2*X'*X;   
%   fprintf('Please install addv and addh.\n');
  end;
 end;
 
catch
  dist=zeros(N);
  tic;
  fprintf('Not enough Memory\nComputing distance line by line ...\n');
  for i=1:N-1
    j=(i+1):N;
    dist(i,j)= ones(1,length(j)).*sum(X(:,i).^2)+sum(X(:,j).^2)-2.*(X(:,j)'*X(:,i))';
    dist(j,i)=dist(i,j)';    
  end;
  fprintf('\n');
end;


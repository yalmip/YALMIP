function [Y,details]=sde(DD,K,varargin)
% [Y,details]=sde(DD,K,pars)
%
%
% DD    SQUARED distance matrix of the input vectors (e.g. euclidean distances)
%
% Optional:
%
% K     number of neighbors 
%
% pars  Parameters
%
%       pars.solver    chooses the SDE solver:
%       pars.solver=0  CSDP   (default)
%       pars.solver=1  SeDuMi
%       pars.solver=2  SDPT3  (very fast - not fully tested)
%
%       pars.slack     switches slack variables on or off
%       pars.slack=0    (default)
%       pars.slack=1   allows distances only to grow 
%       pars.slack=2   allows both growing and shrinking of distances     
%       pars.slack=3   allows distances only to shrink
%                      (in combination with pars.factor=0 leads to
%                      inequality)
%
%       pars.factor    Weight of penalty of Slack variables (only if pars.slac>0)
%                      if pars.factor=0 the equality constraints become
%                      inequalities (direction dependend on pars.slack)                      
%                      default: pars.factor=0.9999 
%
%       pars.inequ=0   default
%       pars.inequ=1   equivalent to pars.slack=3; pars.factor=0;
%                      allos all distances to shrink by setting distance
%                      preserving constraints to inequalities
%
%       pars.verify    Automatically increases K if graph is not connected
%                      (Note: If no K is specified pars.verify=1)
%       pars.verify=1  (default) on
%       pars.verify=0  off
%
%       pars.angles=1  preserves angles in local neighborhood (*DEFAULT*)
%       pars.angles=0  does not preserve angles in local neighborhood 
%
%       pars.repell    Mx2 list of two points that should repell each other
%                      pars.repell=[] (*default*)
%
%       pars.rweight   factor between 0 and 1, of how much the repelling 
%                      should weight compared to the normal objective
%                      function   
%                      pars.rweight=0.3 (*default*)
%
%       pars.kernel    (independent of value) adds the original kernel matrix to the output
%                      (details.K)
%
% Output:
%
% Y         Resulting low dimensional vectors
%           (Rows are dimensions i.e. Y(1:d,:) is d-dimensional embedding 
%
% details   Includes all important details
%           k           size of neighborhood
%           D           all eigenvalues (ascending)
%           info        SDP solver output
%           pars        parameters
%           dual        solution to dual problem
%           inder       order of distance constraints (only interesting for
%                       analysis of slack variables)
%           slack       value of all slack variables
% 
%
% (TIP: From version 1.3 on, it is also possible to pass parameters
% directly as strings
% e.g. [Y,Det]=sde(D,3,'maxiter',20,'solver',0,'inequ'); 
% is equivalent to
%      pars.maxiter=20
%      pars.solver=0
%      pars.inequ=1;  % (if the value is not specified, it is set to 1)
%      [Y,Det]=sde(D,3,pars);
%
% NOTE: sde requires a functionally copy of SeDuMi or CSDP in the path
% (Download CSDP from http://www.nmt.edu/~borchers/csdp.html )
% (and SeDuMi from http://fewcal.kub.nl/sturm/software/sedumi.html )
%
% Example:
%  tt=(linspace(0,1,100).^0.65).*3*pi+pi;  
%  X=[tt.*cos(tt); tt.*sin(tt)];             % generates spiral input data
%  figure;
%  h1=scatter(X(1,:),X(2,:),140,tt,'filled');   % plots input data
%  title('INPUT DATA');
%  set(h1,'MarkerEdgeColor',[0.5 0.5 0.5]);   % draw edges around dots
%  drawnow;
%  axis equal;
%  [Y Det]=sde(distance(X),3,'maxiter',20,'inequ');            % runs SDE (with inequalities)
%  figure;
%  h2=scatter(Y(1,:),Y(2,:),140,tt,'filled');   % plots output data
%  set(h2,'MarkerEdgeColor',[0.5 0.5 0.5]);   % draw edges around dots
%  title('OUTPUT DATA');
%  axis equal; 
%
%
% You can execute this demo with 
% sde('demo');
%
% (Version 1.3)
% (copyright 2004 by Kilian Q. Weinberger:
% http://www.seas.upenn.edu/~kilianw )

if(nargin==0)
  help sde;
  s=input('\n\nShall I run the SDE demo?','s');
  if(length(s)>0 & isequal(upper(s(1)),'Y')) sde('demo');end;
  return;
end;

if(isequal(DD,'demo'))
 tt=(linspace(0,1,100).^0.65).*3*pi+pi;  
 X=[tt.*cos(tt); tt.*sin(tt)];             % generates spiral input data
 figure;
 h1=scatter(X(1,:),X(2,:),140,tt,'filled');   % plots input data
 title('INPUT DATA');
 set(h1,'MarkerEdgeColor',[0.5 0.5 0.5]);   % draw edges around dots
 drawnow;
 axis equal;
 [Y Det]=sde(distance(X),3,'maxiter',20,'inequ');            % runs SDE (with inequalities)
 figure;
 h2=scatter(Y(1,:),Y(2,:),140,tt,'filled');   % plots output data
 set(h2,'MarkerEdgeColor',[0.5 0.5 0.5]);   % draw edges around dots
 title('OUTPUT DATA');
 axis equal;   
 return;
end;



% extract variables
if(length(varargin)>=1)
 if(~isstr(varargin{1}))
    pars=varargin{1};
    for j=1:length(varargin)-1
        varargin{j}=varargin{j+1};
    end;
 end;

 for i=1:nargin-2
  if(isstr(varargin{i}))
    if(i+1>nargin-2 | isstr(varargin{i+1})) val=1;else val=varargin{i+1};end;
    eval(['pars.' varargin{i} '=' sprintf('%i',val) ';']);
  end;
 end;
end;


[t1 t2]=size(DD);
if(t1~=t2) 
    % little check if distance matrix is really square
    error('First argument must be square(!!) distant matrix');
end;

% Set standard parameters
if(~exist('pars') | ~isfield(pars,'factor')) pars.factor=0.9999;end; 
if(~isfield(pars,'slack')) pars.slack=0;end;
if(~isfield(pars,'verify')) pars.verify=1;end;
if(~isfield(pars,'solver')) pars.solver=0;end;
if(~isfield(pars,'angles')) pars.angles=1;end;
if(~isfield(pars,'repell')) pars.repell=[];end;
if(~isfield(pars,'rweight')) pars.rweight=0.3;end;
if(~isfield(pars,'init')) pars.init=0;end;
if(~isfield(pars,'maxiter')) pars.maxiter=50; end;
if(~isfield(pars,'maxiter')) pars.maxiter=50; end;
if(~isfield(pars,'inequ')) pars.inequ=0; end;



% Calculate neighborhoods
if(nargin<2)
  K=neighborsUDD(DD,3);  
else
  if(pars.verify) K=neighborsUDD(DD,K);end;
end;

if(pars.inequ==1)
    pars.slack=3;
    pars.factor=0;
end;


N=length(DD);
[sorted,index] = sort(DD);
neighbors = index(2:(1+K),:);

if(isfield(pars,'ne')) 
  neighbors=pars.ne; % not documented feature to pass as specific
                     % neighborhood matrix as parameter
end;


% Constructing A and b  - the constraints for the SDP
nck=nchoosek(K+1,2);
depth=N*nck*4;

% Initialize matrix A
A=zeros(depth,3);
A(1:end,1)=vec(repmat(1:N*nck,4,1));
b=zeros(N*nck,1);

% Initialize Matrix to check for redundant constraints
DONE=sparse([],[],[],N,N,0);
VALID=ones(N*nck,1);

% just for detailed output, a mapping of input-points > constraints
inder=zeros(depth,2);

% KK = sdpvar(N,N);
% F = set(KK>0) + set(sum(sum(KK))==0);
% for i=1:N
%    ne = neighbors(:,i);
%    pairs=nchoosek([ne;i],2);
%    js=pairs(:,1);
%    ks=pairs(:,2);
%    for ii = 1:length(js)
%        F = F + set(KK(js(ii),js(ii))-2*KK(js(ii),ks(ii))+KK(ks(ii),ks(ii)) == DD(js(ii),ks(ii)));
%    end
% end
% solvesdp(F,-trace(KK));
% KK = double(KK)

% Compute indicies


if 1
allpairs = [];
for i=1:N
   ne = neighbors(:,i);
   pairs=nchoosek([ne;i],2);
   allpairs = [allpairs;pairs];      
end;
allpairs = unique(sort(allpairs,2),'rows');
v1=sub2ind([N N],allpairs(:,2),allpairs(:,2));
v3=sub2ind([N N],allpairs(:,2),allpairs(:,1));
v4=sub2ind([N N],allpairs(:,1),allpairs(:,1)); 
% Setup problem
X = sdpvar(N,N);
s = sdpvar(length(allpairs),1);
if(pars.slack==1)
    F = set(X >= 0) + set(s >= 0);
    F = F + set(sum(sum(X)) == 0) + set(X(v4) - 2*X(v3) + X(v1) + s == DD(v3));
    obj = sum(s)*pars.factor-(1-pars.factor)*trace(X);
else
   F = set(X >= 0);
   F = F + set(sum(sum(X)) == 0) + set(X(v4) - 2*X(v3) + X(v1) == DD(v3));
   obj = -(1-pars.factor)*trace(X);
end

% Interpret in Primal cone
[F,obj] = dualize(F,obj);

% Solve
%solvesdp(F,-obj,sdpsettings('solver','csdp','csdp.maxiter',50));
%X = double(X);
end
% 
Y=[];
details = [];
return

for i=1:N
   ne = neighbors(:,i);
   pairs=nchoosek([ne;i],2);
   
   js=pairs(:,1);
   ks=pairs(:,2);
   
   %find positions in output gram matrix 
   v1=sub2ind([N N],ks,ks);
   v2=sub2ind([N N],js,ks);
   v3=sub2ind([N N],ks,js);
   v4=sub2ind([N N],js,js);

   
   pos=(i-1)*4*nck+1;

   vv=vec([v1 v2 v3 v4]');
   % Eliminate doubles
   mask=DONE(v2)==0;

   if(pars.angles==0)
       mask=mask.*((js==i)+(ks==i));
   end;
   
   VALID((i-1)*nck+1:(i-1)*nck+nck)=mask;
   DONE([v2;v3])=DONE([v2;v3])+1;
   A(pos:pos+4*nck-1,2)=vv;
   A(pos  :4:pos+nck*4-4,3)=A(pos  :4:pos+nck*4-4,3)+1;
   A(pos+1:4:pos+nck*4-3,3)=A(pos+1:4:pos+nck*4-3,3)-1;
   A(pos+2:4:pos+nck*4-2,3)=A(pos+2:4:pos+nck*4-2,3)-1;
   A(pos+3:4:pos+nck*4-1,3)=A(pos+3:4:pos+nck*4-1,3)+1;
            
   
   b(((i-1)*nck+1):(i*nck))= DD(sub2ind([N N],js,ks))';

   inder((i-1)*nck+1:(i*nck),:)=[js ks];
end;
A=sparse(A(:,1),A(:,2),A(:,3));
i=find(VALID>0);
A=A(i,:);
b=b(i,:);
inder=inder(i,:);
fprintf('Redundancy:%d%%\n',round(100-full(length(i))/(length(VALID))*100));


%% Add constraint to center points
A=[ones(1,N^2);A];
b=[0;b];
i=[1;i];



if(isfield(pars,'c')) 
  fprintf('Hidden Feature: Objective Function forced\n');
  c=pars.c;%not documented feature to use weighted 
              % objective function
 else
  c=0-vec(eye(N));
end;

fprintf('Size of A: %ix%i \n\n\n',size(A));

flags.s=N;
flags.l=0;

%Y=[];
%details = [];
%return

%% Add repelling forces
if(length(pars.repell>0))
    [rp,temp]=size(pars.repell);
    c2=zeros(length(c),1);
    js=pars.repell(:,1);
    ks=pars.repell(:,2);   
    v1=sub2ind([N N],ks,ks);
    v2=sub2ind([N N],js,ks);
    v3=sub2ind([N N],ks,js);
    v4=sub2ind([N N],js,js);

    for counter=1:rp
        c2(v1(counter))=1;
        c2(v2(counter))=-1;%
        c2(v3(counter))=-1;%
        c2(v4(counter))=1;
    end;
    
    c=c.*(1-pars.rweight)-c2.*pars.rweight;
end;


if(pars.slack==3)
  % Add Additive Slack Variables
  flags.l=length(i)-1;  
  c=[ones(flags.l,1).*pars.factor;c.*(1-pars.factor)];
  A=[[zeros(1,flags.l); speye(flags.l)] A];
end;


if(pars.slack==2)
  % Add Additive Slack Variables (both additive and sub
  flags.l=(length(i)-1)*2;  
  c=[ones(flags.l,1).*pars.factor;c.*(1-pars.factor)];
  A=[[zeros(1,length(i)-1); speye(length(i)-1)] [zeros(1,length(i)-1); -speye(length(i)-1)] A];
end;

if(pars.slack==1)
  % Add Additive Slack Variables
  flags.l=length(i)-1;  
  c=[ones(flags.l,1).*pars.factor;c.*(1-pars.factor)];
  A=[[zeros(1,flags.l); -speye(flags.l)] A];
end;

if(pars.slack & 1==0)
  % Add Multiplicative Slack Variables
  flags.l=length(i)-1;  
  c=[ones(flags.l,1).*pars.factor;c.*(1-pars.factor)];
  A=[[zeros(1,flags.l); -sparse(spdiags(b(2:end),0,length(b)-1,length(b)-1))] A];
end;


Kneigh=K;

% Y=[];
% details=[];
% return

switch(pars.solver)
    case{0},
      OPTIONS.maxiter=50;%pars.maxiter;
      [x d z  info]=csdp(A,b,c,flags,OPTIONS);
      K=mat(x(flags.l+1:flags.l+flags.s^2));
      slack=x(1:flags.l);
    case{1},
       OPTIONS.maxiter=pars.maxiter;
       [k d info]=sedumi(A,b,c,flags,OPTIONS);
      K=mat(k(end-flags.s^2+1:end));
      slack=k(1:end-flags.s^2);
      %keyboard;
 
 case{2},
       fprintf('converting into SQLP...');
      [blk,A,c,b]=convert_sedumi(A',b,c,flags);
      fprintf('DONE\n');
       OPTIONS.vers=1;
       OPTIONS.gam=0;
       OPTIONS.predcorr=1;
       OPTIONS.expon=[1 1];
       OPTIONS.gaptol=1.0000e-08;
       OPTIONS.inftol=1.0000e-08;
       OPTIONS.steptol=1.0000e-06;
       OPTIONS.maxit=pars.maxiter;
       OPTIONS.printyes=1;
       OPTIONS.scale_data=0;
       OPTIONS.randnstate=0;
       OPTIONS.spdensity=0.5000;
       OPTIONS.rmdepconstr=0;
       OPTIONS.cachesize=256;
       OPTIONS.smallblkdim=15;
    
      if(pars.init==1)  
         fprintf('Initializing ...');
          % Initialize with input points
        [X0,y0,Z0]=infeaspt(blk,A',c,b); 
        X0{length(X0)}=getgram(DD)+eye(N).*0.001;
        fprintf('DONE\n');
        [obj,Kt,d,z]=sqlp(blk,A,c,b,{},X0,y0,Z0);
      else
        OPTIONS.maxit=pars.maxiter;
        [obj,Kt,d,z]=sqlp(blk,A,c,b,OPTIONS);
    end;
      K=Kt{length(Kt)};
      info=0;
      if(pars.slack)slack=Kt{1}(:);end;
end;
%keyboard;
[V D]=eig(K);

V=V*sqrt(D);
D=diag(D);

Y=(V(:,end:-1:1))';

if(isfield(pars,'kernel')) 
 details.K=K;
end;
details.k=Kneigh;
details.D=D;
details.info=info;
details.pars=pars;
details.dual=d;
details.inder=inder;


if(pars.slack)details.slack=slack;end;

return;












 






%%% SIMPLE BUT USEFUL UTIL FUNCTIONS

function neighbors=getneighborsUDD(DD,K);
ne=getneighborsDD(DD,K);
N=length(DD);

for i=1:N
    neighbors{i}=[];
end;

for i=1:N
 neighbors{i}=merge(sort(neighbors{i}),sort(ne(:,i)));
 for j=1:K
    neighbors{ne(j,i)}=merge(neighbors{ne(j,i)}, i);
 end;
end;



function result=merge(x,y);
result=unique([x(:);y(:)]);


function v=vec(M);
v=M(:);


function k=neighborsUDD(DD,K);
N=length(DD);
if(nargin<2)
    K=2;
end;
k=K;
while(k<N & (1-connectedUDD(DD,k)))
    k=k+1;
    fprintf('Trying K=%i \n',k);
end;



function result=connectedUDD(DD,K);
% result = connecteddfs (X,K)
%
% X input vector
% K number of neighbors
%
% Returns: result = 1 if connected 0 if not connected

if(nargin<2)
    fprintf('Number of Neighbors not specified!\nSetting K=4\n');
    K=4;    
end;
N=length(DD);


ne=getneighborsUDD(DD,K);

 
maxSize=0;
for i=1:N
    if(length(ne{i})>maxSize) maxSize=length(ne{i});end;
end;
neighbors=ones(maxSize,N);
for i=1:N
    neighbors(1:length(ne{i}),i)=ne{i};    
end;
oldchecked=[];
checked=[1];

while((size(checked)<N) & (length(oldchecked)~=length(checked)))  
 next=neighbors(:,checked);
 next=transpose(sort(next(:)));
 next2=[next(2:end) 0];
 k=find(next-next2~=0);
 next=next(k);

 oldchecked=checked; 
 checked=neighbors(:,next(:));
 checked=transpose(sort(checked(:)));
 checked2=[checked(2:end) 0];
 k=find(checked-checked2~=0);
 checked=checked(k);
end;
result=(length(checked)==N);



function X = mat(x,n)
 if nargin < 2
   n = floor(sqrt(length(x)));
   if (n*n) ~= length(x)
     error('Argument X has to be a square matrix')
   end
 end
 X = reshape(x,n,n);

 
 
 
function neighbors=getneighborsDD(DD,K);
N=length(DD);
[sorted,index] = sort(DD);
neighbors = index(2:(1+K),:);

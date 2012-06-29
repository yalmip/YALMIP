function [Y,details]=sdeIncXL(X,k,FirstSet,d,pars);
%
% [Y,details]=sdeInc(X,k,FirstSet,d,pars);
%
%
% Required parameters:
%
% X         matrix containing the input vectors as columns
% k         size of local neighborhood
% FirstSet  Size of Subsample 
%           (only for speedup purpose 
% d         dimension of output
%
%
% Optional parameters:
%
% pars:
%              pars.slack=1 Use slack for first set (*DEFAULT*)
%              pars.slack=0 Don't use slack for first set
%               
%              pars.solver=0 Use CSDP for non-first sets (*DEFAULT*)
%              pars.solver=1 Use SeDuMi for non-first sets 
%              pars.solver=2  SDPT3  (very fast - not fully tested)
%
%              pars.verify    Automatically increases K if graph is not connected
%                      (Note: If no K is specified pars.verify=1)
%              pars.verify=1  (default) on
%              pars.verify=0  off
%
%
%              pars.factor    Weight of Slack variables (only if pars.slack=1)
%                             the very closer to 1 ther harder do the slack
%                             variables get
%              pars.factor=0.9999 (default)
%
%              pars.angles=1  preserves angles in local neighborhood (*DEFAULT*)
%              pars.angles=0  does not preserve angles in local neighborhood 
%
%              pars.repell     Mx2 list of two points that should repell each other
%                              pars.repell=[] (*default*)
%
%              pars.rweight    factor between 0 and 1, of how much the repelling 
%                              should weight compared to the normal objective
%                              function   
%                              pars.rweight=0.3 (*default*)
%
%              pars.save       Filename to save temporare states
%                              default='tempsave'
%
%              pars.continue   Filename to continue from temporare states
%                              default=''
%
%              pars.i          initial subset (default random)
%
% Output:
% 
% Y            Matrix with output vectors 
% details
% details.SY   Eigenvectors of subsample
% details.SD   Eigenvalues of subsample
% details.time time needed for computation
%
% (version 1.2)
% (copyright 2004 by Kilian Q. Weinberger:
% http://www.seas.upenn.edu/~kilianw )
[D N]=size(X);
i=randperm(N); % bug in matlab

if(~exist('pars') | ~isfield(pars,'factor')) pars.factor=0.9999;end; 
if(~isfield(pars,'slack')) pars.slack=1;end;
if(~isfield(pars,'verify')) pars.verify=1;end;
if(~isfield(pars,'solver')) pars.solver=0;end;
if(~isfield(pars,'repell')) pars.repell=[];end;
if(~isfield(pars,'rweight')) pars.rweight=0.3;end;
if(~isfield(pars,'save')) pars.save='tempsave.mat';end;
if(~isfield(pars,'i')) pars.i=randperm(N);end;

A=sparse([],[],[],d*(d+1)/2,(d+1)^2);
acounter=1;dim=d;b=[];
for i=1:dim
      for j=i:dim
       a=zeros(dim+1);
       a(i,j)=1;
       A(acounter,:)=reshape(a,1,(dim+1)^2);
       acounter=acounter+1;
       b=[b;i==j];
      end;
end;
pars.stampA=A;
pars.stampb=b;
clear('A');


pars


if(~isfield(pars,'continue')) 

fprintf('Dividing problem into subproblems ...');
 % this line is executed twice due to a bug in Matlab 6
i=pars.i;
if(length(pars.repell)>0)
   repeller=unique(pars.repell(:));
%   keyboard;
   [temp,temp2,reppos]=intersect(repeller,i);
%   [temp,reppos]=ismember(repeller,i)
   if(length(reppos>0))
    i=[repeller' setdiff(i,repeller)];
    
    j=zeros(1,length(X));
    j(repeller)=1:length(repeller);
    pars.repell(:)=j(pars.repell(:));
   end; 
end;
firsti=i(1:FirstSet);
fprintf('done\n');

    
fprintf('computing distances  ....');
X2 = sum(X.^2);
DMfirsti = repmat(X2(firsti),FirstSet,1)+repmat((X2(firsti))',1,FirstSet)-2*(X(:,firsti)'*X(:,firsti));
%clear('X2');
fprintf('done\n');


if(pars.verify)
 fprintf('Checking if graph is connected ...');
 k=neighborsUDD(DMfirsti,k);
 fprintf('done\n');
end;

[V,DD]=sde(DMfirsti,k,pars);

Y=V(1:d,:);

index=firsti;clear('firsti');
jj=ones(1,N);
jj(index)=0;
jj=find(jj);


if(~strcmp(pars.save,''))
      fprintf('Saving temporary file...');
      save(pars.save);
      fprintf('done\n');
end;

fprintf('Sorting additional points ...\n');
distances=zeros(1,length(jj));
for tt1=1:length(distances);
    x=X(:,jj(tt1));
    distances(tt1)=min(X2(index)'-((X(:,index))'*x).*2+x'*x);  
    if(mod(tt1,1000)==0)  fprintf('Points processed: %i/%i\n',tt1,length(distances)); end;
end;
[temp,tt1]=sort(distances);
jj=jj(tt1);

fprintf('done\n');


ScaleFactor=max(max(abs(Y)));
Y=Y./ScaleFactor;
%DM=DM./(ScaleFactor^2);

err=[];
shift=zeros(d,1);
precomputed=0;

else
    load(pars.continue);    
    [temp,precomputed]=size(Y);
end;

tic;
while(length(jj)>0)
  
  [temp ntotal]=size(jj);
  
%distances=DM(index,jj);
%[temp iio]=sort(min(distances));s
  ii=1;
  
   x=X(:,jj(ii));


   distances=X2(index)'-((X(:,index))'*x).*2+x'*x;

   
   disttone=zeros(1,k);
   ne=zeros(1,k);
   for te=1:k
    [tt1 tt2]=min(distances);
    disttone(te)=tt1;
    ne(te)=tt2;
    distances(tt2)=inf;    
   end;
   % [disttone ne]=sort(distances(:,ii));

   

   [y er]=sdeAddPoint3(Y(:,ne)-repmat(shift,1,k),disttone./ScaleFactor^2,pars);
   Y=[Y,y+shift];
%   XX=[XX x];
   err=[err er.err]; 
    
%  save('temp4');
 index=[index jj(ii)];

 jj=jj(2:end);
 
% Xle=X(:,jj);

 % centralize Y
 %Y=Y-repmat(y./(N-ntotal+1),1,N-ntotal+1);
 shift=shift+y./(N-ntotal+1);
 lindex=N-ntotal+1;
 fprintf('%i ',lindex);
 
 if(mod(lindex,100)==0)
  timeremaining(length(index)-precomputed,length(jj));
  if(mod(lindex,1000)==0 & ~strcmp(pars.save,''))
      save(pars.save);
  end;
 end;
end;

details.SY=Y(:,1:FirstSet).*ScaleFactor;	
details.index=index;
details.SD=DD;
details.time=toc;
details.err=err;
details.k=k;
Y(:,index)=Y.*ScaleFactor;









function [y,details]=sdeAddPoint3(Anchors,distances,pars);

if(~exist('pars') | ~isfield(pars,'factor')) pars.factor=0.9999;end; 
if(~isfield(pars,'solver')) pars.solver=0;end;
if(~isfield(pars,'rweight')) pars.rweight=0;end;


  
[dim,k]=size(Anchors);

% Construct A
% First fix the identity matrix in Zxs
%a=[ones(dim,1); 1];
%A=sparse([],[],[],dim*(dim+1)/2+k,(dim+1)^2);
b=[];

%acounter=1;

  b=pars.stampb;
 A=pars.stampA;
 [temp1 temp2]=size(A);
 acounter=temp1+1;


% Now add anchor constraints
for i=1:k
  a=[Anchors(:,i); -1];
     
  A(acounter,:)=[reshape(a*a',1,(length(a))^2)];
  acounter=acounter+1;

  b=[b; distances(i)];
end;


[constraints ll]=size(A);
A=[[zeros(constraints-k+1,2*(k-1));[eye(k-1)  -1*eye(k-1)]] A];

K.l=2*(k-1);
K.s=dim+1;

if(pars.rweight==1)
  c=[ones(1,K.l) zeros(1,(K.s)^2)];  
else
  c=[ones(1,K.l)*pars.factor zeros(1,(K.s)^2-1) -1*(1-pars.factor)  ];
end;



%if(pars.solver)
switch(pars.solver)  
 case{0},
  pars2.printlevel=0;
  [kk temp temp info]=csdp(A,b,c,K,pars2);
  M=reshape(kk(K.l+1:K.l+K.s^2),K.s,K.s);
 case{1},
  pars2.fid=0;
  [kk temp info]=sedumi(A,b,c,K,pars2);
  M=reshape(kk(end-K.s^2+1:end),K.s,K.s);
 case{2},
%  keyboard;
  %fprintf('converting into SQLP...');
  [blk,A,c,b]=convert_sedumi(A',b,c,K);
  %fprintf('DONE\n');
  sqlparameters;
  OPTIONS.verbose=0;
  [obj,kkt,temp,temp]=sqlp(blk,A,c,b,OPTIONS);
  M=kkt{length(kkt)};
  %keyboard;
  %M=reshape(kk(K.l+1:K.l+K.s^2),K.s,K.s);
  info=0;
end; 

y=M(1:dim,end);
G=M(end,end);
%fprintf('Trace Error: %d\n',abs(G-y'*y));
details.err=abs(G-y'*y);


















%%% SIMPLE BUT USEFUL UTIL FUNCTIONS


function timeremaining(sofar,left);
  time=toc/(sofar)*left;
  days=floor(time/(60^2*24));
  hours=floor((time-days*60^2*24)/60^2);
  minutes=floor((time-hours*60^2)/60);
  seconds=floor(time-minutes*60);

  fprintf('\nTime remaining: ');
  if(days>0) fprintf('%i Day(s) ',days);end;
  if(hours>0  ) fprintf('%i Hour(s) ',hours);end;
  if(minutes>0 & days==0) fprintf('%i Minute(s) ',minutes);end;
  if(seconds>0 & hours==0 & days==0) fprintf('%i Second(s) ',seconds);end;
  if(time==0) fprintf('You are DONE!!');end;
  fprintf('\n');  


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
%  if(length(oldchecked)==length(checked))
%      prod(double(checked==oldchecked));     
%  end;
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
% PAIRWISE DISTANCES
N=length(DD);
% NEIGHBORS
[sorted,index] = sort(DD);
neighbors = index(2:(1+K),:);

    








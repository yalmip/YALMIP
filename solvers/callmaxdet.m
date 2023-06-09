function output = callmaxdet(interfacedata)

[F_struc,F_blksz,G_struc,G_blksz] = sedumi2maxdet(interfacedata.F_struc,interfacedata.K); 

c = interfacedata.c;
options = interfacedata.options;
x0 = interfacedata.x0;
onlyfeasible = (nnz(c)==0) & isempty(G_blksz);

if isempty(G_blksz)
    G_blksz = 1;
    G_struc = [1 zeros(1,length(c))];
else
end

if isempty(F_blksz)
    x = zeros(length(c),1)+NaN;
    problem = 11;
    infostr = yalmiperror(problems,'MAXDET failed since there are no inequality constraints in F(x)');
    solveroutput = [];
    solvertime = 0;
    D_struc = [];  
end

problem = 0;
solvertimephase1=0;
D_struc = [];
if isempty(x0)
    solvertimephase1 = tic;
    showprogress('Calling MAXDET/phase1',options.showprogress);
    [x0,z0,w0,problem,infostr,solveroutput] = callmaxdetphase1(full(F_struc),full(F_blksz), full(G_struc),full(G_blksz), full(c), options);
    solvertimephase1 = toc(solvertimephase1);
end
if (problem~=0) | (onlyfeasible==1)
    if isempty(x0)
        x = repmat(nan,length(c),1);
    else
        x = x0;	
    end
    solvertime = solvertimephase1;
else
    z0=zeros(size(F_struc,1),1);
    w0=zeros(size(G_struc,1),1);	
    solvertime = tic;
    showprogress('Calling MAXDET',options.showprogress);
    if options.verbose==0
        evalc('[x,Z,W,ul,hist,infostr]=maxdet(full(F_struc),F_blksz,full(G_struc), G_blksz,c'',x0,z0,w0,options.maxdet.AbsTol,options.maxdet.RelTol,options.maxdet.gam,options.maxdet.NTiters);');
    else
        [x,Z,W,ul,hist,infostr]=maxdet(full(F_struc),F_blksz,full(G_struc), G_blksz,c',x0,z0,w0,options.maxdet.AbsTol,options.maxdet.RelTol,options.maxdet.gam,options.maxdet.NTiters);
    end
    solvertime = toc(solvertime)+solvertimephase1;
    
    D_struc = [Z;W];
    
    if options.savesolveroutput
        solveroutput.x = x;
        solveroutput.Z = Z;
        solveroutput.W = W;
        solveroutput.ul = ul;
        solveroutput.hist = hist;
        solveroutput.infostr = infostr;
    else
        solveroutput = [];
    end
    
    switch  infostr
        case {'target reached','relative accuracy reached','absolute accuracy reached','absolute tolerance reached','relative tolerance reached'}
            problem = 0;
        case 'maximum Newton iteration exceeded'
            problem = 3;
        otherwise
            problem = 1;
    end    
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,[],solveroutput,solvertime);

function [x0,z0,w0,problem,infostr,solveroutput] = callmaxdetphase1(F_struc,F_blksz, G_struc,G_blksz, c, options);

try
    if options.verbose==0
        evalc('[x0,z0,w0,ul] = phase1(full(F_struc),F_blksz,full(G_struc), G_blksz,options.maxdet.gam,options.maxdet.AbsTol,options.maxdet.RelTol,options.maxdet.NTiters);');
    else
        [x0,z0,w0,ul] = phase1(full(F_struc),F_blksz,full(G_struc), G_blksz,options.maxdet.gam,options.maxdet.AbsTol,options.maxdet.RelTol,options.maxdet.NTiters);
    end
    if ul(1)<0
        problem = 0;
    elseif ul(2)>=0
        problem = 1;
    else
        problem = 8;
    end
    %  % Problems currently detected outside
%  problem = (ul(1)>0);

% In case of problems, phase1 outputs the extended vector instead
% of last primal iterate
if problem
    x0 = x0(1:length(c));
end

catch
    problem = 9;
    x0 = [];
    z0 = [];
    w0 = [];
    ul = [];
end
infostr = yalmiperror(problem,'MAXDET/phase1');
if options.savesolveroutput
    solveroutput.x = x0;
    solveroutput.Z = z0;
    solveroutput.W = w0;
    solveroutput.ul = ul;
    solveroutput.infostr = infostr;
else
    solveroutput = [];
end

% CODE TAKEN DIRECTLY FROM MAXDET/PHASE1
% THE REASON FOR HAVING A COPY HERE IS THAT 
% A FILE NAMED PHASE1 EXIST FOR SP ALSO
%MAXDET, version ALPHA, April 1996.

%COPYRIGHT 1996 SHAO-PO WU, LIEVEN VANDENBERGHE AND STEPHEN BOYD
%Permission to use, copy, modify, and distribute this software for 
%any purpose without fee is hereby granted, provided that this entire 
%notice is included in all copies of any software which is or includes
%a copy or modification of this software and in all copies of the 
%supporting documentation for such software.
%This software is being provided "as is", without any express or 
%implied warranty.  In particular, the authors do not make any
%representation or warranty of any kind concerning the merchantability
%of this software or its fitness for any particular purpose.
function [x,Z,W,ul] = phase1(F,F_blkszs,G,G_blkszs,gam,abstol,reltol,NTiters);

m = size(F,2)-1;
if (m ~= size(G,2)-1)
    error('F and G must have the same number of columns.');
end
if (size(F,1) ~= sum(F_blkszs.^2)) 
    error('Dimensions of F do not match F_blkszs.');
end;
if (size(G,1) ~= sum(G_blkszs.^2)) 
    error('Dimensions of G do not match G_blkszs.');
end;

% mineigF is the smallest eigenvalue of F_0
mineigF = 0.0;
k=0; for n=F_blkszs,
    mineigF = min(mineigF, min(eig(reshape(F(k+[1:n*n],1),n,n))));
    k=k+n*n;   % k = sum n_i*n_i 
end;
% mineigG is the smallest eigenvalue of G_0
mineigG = 0.0;
k=0; for n=G_blkszs,
    mineigG = min(mineigG, min(eig(reshape(G(k+[1:n*n],1),n,n))));
    k=k+n*n;   % k = sum n_i*n_i 
end;

% eyeF is the identity
eyeF = zeros(size(F,1),1);  
k=0; for n=F_blkszs,
    eyeF(k+[1:n*n]) = reshape(eye(n),n*n,1);   % identity
    k=k+n*n;   % k = sum n_i*n_i 
end;
% eyeG is the identity
eyeG = zeros(size(G,1),1);  
k=0; for n=G_blkszs,
    eyeG(k+[1:n*n]) = reshape(eye(n),n*n,1);   % identity
    k=k+n*n;   % k = sum n_i*n_i 
end;

% initial x0 
x0 = [zeros(m,1); max(-1.1*min(mineigF,mineigG), 1)];

% linear objective
c = [zeros(m,1); 1];

% call maxdet
[x,Z,W,ul,hist,infostr]=maxdet([F,eyeF; G,eyeG],...
    [F_blkszs(:)',G_blkszs(:)'],...
    [1 zeros(1,m+1)],1,c,x0,...
    zeros(size(F,1)+size(G,1),1),0,...
    abstol,reltol,gam,NTiters);

% prepare output
x = x(1:m);
W = Z(size(F,1)+1:size(F,1)+size(G,1));
Z = Z(1:size(F,1));
%if ul(1)<0
%  infostr = 'feasible';
%elseif ul(2)>=0
%  infostr = 'infeasible';
%else
%  infostr = 'feasibility cannot be determined';
%end

function [F_struc,blksz] = deblock(F_struc,blksz);
X = any(F_struc(end-blksz^2+1:end,:),2);
X = reshape(X,blksz,blksz);
[v,dummy,r,dummy2]=dmperm(X);
blks = diff(r);

logt = F_struc;

newlogt = [];
for i = 1:size(logt,2)
    temp = reshape(logt(:,i),blksz,blksz);
    temp = temp(v,v);
    newlogt = [newlogt temp(:)];
end
logt = newlogt;

pattern = [];
for i = 1:length(blks)
    pattern = blkdiag(pattern,ones(blks(i)));
end

F_struc = logt(find(pattern),:);
blksz = blks;
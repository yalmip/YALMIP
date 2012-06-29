%SEDUMI2SDPT3 Internal function, obsolete

%%*******************************************************************
%%  Converts from SeDuMi format.
%%
%%  [blk,A,C,b] = sedumi2sdpt3(c,A,b,K)
%%
%%  Input: fname.mat = name of the file containing SDP data in
%%                     SeDuMi format.
%%
%% Important note: Sedumi's notation for free variables "K.f"
%%                 is coded in SDPT3 as blk{p,1} = 'u', where
%%                 "u" is used for unrestricted variables.
%%
%% Ripped from:
%% SDPT3: version 3.0
%% Copyright (c) 2000 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************************

function [blk,Avec,C,b,oldKs] = sedumi2sdpt3(c,A,b,K,smallblkdim)

if (size(c,1) == 1), c = c'; end;
if (size(b,1) == 1), b = b'; end;
if (norm(A,'fro') > 0) & (size(A,2) == length(b)); At = A; end
if (norm(At,'fro')==0), At = A'; end;
[nn,mm] = size(At); if (max(size(c)) == 1); c = c*ones(nn,1); end;
if ~isfield(K,'f'); K.f = 0; end
if ~isfield(K,'l'); K.l = 0; end
if ~isfield(K,'q'); K.q = 0; end
if ~isfield(K,'s'); K.s = 0; end
if (K.f == 0) | isempty(K.f); K.f = 0; end;
if (K.l == 0) | isempty(K.l); K.l = 0; end;
if (sum(K.q) == 0) | isempty(K.q); K.q = 0; end
if (sum(K.s) == 0) | isempty(K.s); K.s = 0; end
%%
%%
%%
m = length(b);
rowidx = 0;  idxblk = 0;
if ~(K.f == 0)
    len = K.f;
    idxblk = idxblk + 1;
    blk{idxblk,1} = 'u'; blk{idxblk,2} = K.f;
    Avec{idxblk,1} = At(rowidx+[1:len],:);
    C{idxblk,1} = c(rowidx+[1:len]);
    rowidx = rowidx + len;
end
if ~(K.l == 0)
    len = K.l;
    idxblk = idxblk + 1;
    blk{idxblk,1} = 'l'; blk{idxblk,2} = K.l;
    Avec{idxblk,1} = At(rowidx+[1:len],:);
    C{idxblk,1} = c(rowidx+[1:len]);
    rowidx = rowidx + len;
end
if ~(K.q == 0)
    len = sum(K.q);
    idxblk = idxblk + 1;
    blk{idxblk,1} = 'q'; blk{idxblk,2} = K.q;
    Avec{idxblk,1} = At(rowidx+[1:len],:);
    C{idxblk,1} = c(rowidx+[1:len]);
    rowidx = rowidx + len;
end
oldKs = [];
if ~(K.s == 0)
    % Avoid extracting rows!
    At = At';
    smblkdim = smallblkdim;
    blksize = K.s;
    if size(blksize,2) == 1; blksize = blksize'; end
    blknnz = [0 cumsum(blksize.*blksize)];
    deblkidx = find(blksize > smblkdim);
    if ~isempty(deblkidx)
        for p = 1:length(deblkidx)
            idxblk = idxblk + 1;
            n = blksize(deblkidx(p));
            oldKs = [oldKs deblkidx(p)];
            pblk{1,1} = 's'; pblk{1,2} = n;
            blk(idxblk,:) = pblk;
            Atmp = At(:,rowidx+blknnz(deblkidx(p))+[1:n*n])';
            Avec{idxblk,1} = sparse(n*(n+1)/2,m);
            
            warning_yes = 1;
            if 1
                Avec{idxblk,1} = yalmipsvec(Atmp,n);
            else                
                for k = 1:m
                    Ak = mexmat(pblk,Atmp(:,k),1);
                    Avec{idxblk,1}(:,k) = svec(pblk,Ak,1);
                end
            end            
            Ctmp = c(rowidx+blknnz(deblkidx(p))+[1:n*n]);
            Ctmp = mexmat(pblk,Ctmp,1);
            C{idxblk,1} = 0.5*(Ctmp+Ctmp');
        end
    end
    spblkidx = find(blksize <= smblkdim);
    if ~isempty(spblkidx)
        pos = []; len = 0;
        for p = 1:length(spblkidx)
            n = blksize(spblkidx(p));            
            oldKs = [oldKs spblkidx(p)];
            len = len + n*(n+1)/2;
            pos = [pos, rowidx+blknnz(spblkidx(p))+[1:n*n]];
        end
        idxblk = idxblk + 1;
        blk{idxblk,1} = 's';  blk{idxblk,2} = blksize(spblkidx);
        Avec{idxblk,1} = sparse(len,m);      
        Atmp = At(:,pos)';
        warning_yes = 1;
        for k = 1:m
            Ak = mexmat(blk(idxblk,:),sparse(full(Atmp(:,k))),1);
            Avec{idxblk,1}(:,k) = svec(blk(idxblk,:),Ak,1);
        end
        Ctmp = c(pos);
        Ctmp = mexmat(blk(idxblk,:),Ctmp,1);
        C{idxblk,1} = 0.5*(Ctmp+Ctmp');
    end
end

function x = yalmipsvec(X,n)

Y = reshape(1:n^2,n,n)';
d = diag(Y);
Y = tril(Y);
Y = (Y+Y')-diag(sparse(diag(Y)));
[uu,oo,pp] = unique(Y(:));
X = X*sqrt(2);
X(d,:)=X(d,:)/sqrt(2);
x = X(uu,:);



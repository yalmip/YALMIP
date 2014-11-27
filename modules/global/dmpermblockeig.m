function [V,D,permutation,failure] = dmpermblockeig(X,switchtosparse)
    
[permutation,aux1,aux2,blocks] = dmperm(X+speye(length(X)));
Xpermuted = X(permutation,permutation);

V = [];
D = [];
V = zeros(size(X,1),1);
top = 1;
left = 1;
anycholfail = 0;
failure = 0;

for i = 1:length(blocks)-1
    Xi = Xpermuted(blocks(i):blocks(i+1)-1,blocks(i):blocks(i+1)-1);
    [R,fail] = chol(Xi);
    anycholfail = anycholfail | fail;
    if fail
        if length(Xi) >= switchtosparse           
            [vi,di,eigfail] = eigs(Xi,5,'SA');            
            if eigfail || isempty(di)
                res = 0;
                for j = 1:size(vi,2)
                    res(j) = norm(Xi*vi(:,j)-vi(:,j)*di(j,j));
                end
                % We only trust these
                notfailed = abs(res) <= 1e-12;
                vi = vi(:,notfailed);
                di = di(notfailed,notfailed);
                if length(vi) == 0
                    [vi,di,eigfail] = eigs(sparse(Xi),25,'SA');                    
                    if eigfail
                        res = 0;
                        for j = 1:size(vi,2)
                            res(j) = norm(Xi*vi(:,j)-vi(:,j)*di(j,j));
                        end
                        % We only trust these
                        notfailed = abs(res) <= 1e-12;
                        vi = vi(:,notfailed);
                        di = di(notfailed,notfailed);
                    end
                end
            end
        else
            [vi,di] = eig(full(Xi));
        end
        for j = 1:length(di)
            if di(j,j)<=0
                V(top:top+length(Xi)-1,left)=vi(:,j);
                left = left + 1;
                D = blkdiag(D,di(1,1));
            end
        end
    end
    top = top + length(Xi);
end


if (anycholfail && isempty(V)) || (anycholfail && all(diag(D)>0)) 
    % OK, we have a problem. The Cholesky factorization failed for some of
    % the matrices, but yet no eigenvalue decomposition revealed a negative
    % eigenvalue (due to convergence issues in the sparse eigs)
    failure = 1;
end

function [vi,di,eigfail] = eigband(X,m)

eigfail = 0;
if nargin == 1
    m = 5;
end

if m > 0
    r = symrcm(X);
    Z = X(r,r);
    n = length(Z);
    bw = yalmipbandwidth(Z);
   % spy(Z);drawnow
    if bw > n/3 || n < 200
        [vi,di] = eig(full(X));
        return
    end
    mid = ceil(n/2);
    Z1 = Z(1:mid+bw,1:mid+bw);
    Z2 = Z(mid-bw:end,mid-bw:end);
    [v1,d1,eigfail] = eigband(Z1,m-1);
    [v2,d2,eigfail] = eigband(Z2,m-1);
    vi = blkdiag(v1,v2);
    di = blkdiag(d1,d2);
    i1 = find(diag(d1)<0);
    i2 = find(diag(d2)<0);
    vi = zeros(n,0);
    for i = 1:length(i1)
        vi(1:n/2+bw,end+1) = v1(:,i1(i));
    end
    for i = 1:length(i2)
        vi(n/2-bw:end,end+1) = v2(:,i2(i));
    end
    di = blkdiag(d1(i1,i1),d2(i2,i2));
    [~,ir] = ismember(1:length(r),r);
    vi = vi(ir,:);
    return
end

r = symrcm(X);
Z = X(r,r);
n = length(Z);
bw = yalmipbandwidth(Z);
mid = ceil(n/2);
Z1 = Z(1:mid+2*bw,1:mid+2*bw);
Z2 = Z(mid-2*bw:end,mid-2*bw:end);
[v1,d1] = eig(full(Z1));
[v2,d2] = eig(full(Z2));
i1 = find(diag(d1)<0);
i2 = find(diag(d2)<0);
vi = zeros(n,0);
for i = 1:length(i1)
    vi(1:n/2+bw,end+1) = v1(:,i1(i));
end
for i = 1:length(i2)
    vi(n/2-bw:end,end+1) = v2(:,i2(i));
end
di = blkdiag(d1(i1,i1),d2(i2,i2));
[~,ir] = ismember(1:length(r),r);
vi = vi(ir,:);




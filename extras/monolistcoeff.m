function m = monolistcoeff(n,d,summax)
%MONOLISTCOEFF  Internal function used in SOS programs

if length(d)==1
    d = repmat(full(d),1,full(n));
end

if nargin==2
    summax = (sum(d));
end

if 0
    dmax = sum(d);%max(full(d));%sum(d);
    base = eye(n);
    m = base;
    for i = 1:dmax-1
        temp=[];
        for k = 1:n
            temp = [temp;m+repmat(base(k,:),size(m,1),1)];
        end
        ii=find(~any(temp-repmat(d,size(temp,1),1)>0,2));
        temp=temp(ii,:);
        ii=find(~any(sum(temp,2)-summax>0,2));
        temp=temp(ii,:);

        m = [m;temp];
        [ii,jj,kk]=uniquesafe(m,'rows');
        m = m(jj,:);
    end
    % Add constant
    m = [zeros(1,n);m];
else
    m = (0:d(1))';
    for i = 2:n
        z = (0:d(i))';
        m = [kron(m,ones(size(z,1),1)) kron(ones(size(m,1),1),z)];
        [ii,jj,kk]=uniquesafe(m,'rows');
        m = m(jj,:);
        m = m(sum(m,2)<=summax,:);
    end    
end
m = full(m);
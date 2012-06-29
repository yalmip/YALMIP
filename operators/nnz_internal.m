function F = nnz_internal(z,x,direction)

switch direction
    case 0
        [nx,mx] = size(x);
        if issymmetric(x)
            d = binvar(nx,nx);  % d == 1 means x(i,j) can be non-zero
        else
            d = binvar(nx,mx,'full');
        end
        x = reshape(x,nx*mx,1);
        d = reshape(d,nx*mx,1);
        [M,m] = derivebounds(x);
        F = set(m.*d <= x <= M.*d) + set(z == sum(sum(d)));
        F = F + set(0 <=z <= nx*mx);

    case 1
        [nx,mx] = size(x);
        if issymmetric(x)
            du = binvar(nx,nx);  % du == 1 means x(i,j) > 0
            dd = binvar(nx,nx);  % dd == 1 means x(i,j) < 0
        else
            du = binvar(nx,mx,'full');
            dd = binvar(nx,mx,'full');
        end
        x = reshape(x,nx*mx,1);
        du = reshape(du,nx*mx,1);
        dd = reshape(dd,nx*mx,1);
        [M,m] = derivebounds(x);
        fixedzeros = find(M==m & m==0);       
        positive = find(m>0);
        negative = find(M<0);
        left = setdiff(1:nx*mx,[fixedzeros;positive;negative]);
        
        m = m(left);
        M = M(left);
        x = x(left);
        du = du(left);
        dd = dd(left);
        du(M<=0) = 0;
        dd(m>=0) = 0;
        F = [];
        eps=1e-3;
        F = F + set(sum([du dd],2) <= 1);
        F = F + set(z == length(negative)+length(positive) + ((sum(du) + sum(dd))));
        F = F + set(x >= eps+(m-eps).*(1-du));
        F = F + set(x <= eps+(M-eps).*du);
        F = F + set(x <=-eps+(M+eps).*(1-dd));
        F = F + set(x >= -eps+(m+eps).*dd);
        F = F + set(m.*(du+dd) <= x <= M.*(dd+du));
        F = F + set(0 <= z <=nx*mx);
        
        %F = F + set(x > 0+(m-0).*(1-du));
        % F = F + set(x < -0+(0+M).*(1-dd));
        %F = F + set(0 <= z <=nx*mx);
        
      %  F = set(m.*(du+dd) <= x <= M.*(du+dd)) + set(sum([du dd],2) <= 1);
      %  F = F + set(z == ((sum(du) + sum(dd))));
      %  F = F + set(x > 1e-5+(-1+m).*(1-du));
      %  F = F + set(x < -1e-5+(1+M).*(1-dd));
      %  F = F + set(0 <= z <=nx*mx);
    otherwise
end


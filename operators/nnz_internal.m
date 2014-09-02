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
        F = (m.*d <= x <= M.*d) + (z == sum(sum(d)));
        F = F + (0 <=z <= nx*mx);

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
        F = F + (sum([du dd],2) <= 1);
        F = F + (z == length(negative)+length(positive) + ((sum(du) + sum(dd))));
        F = F + (x >= eps+(m-eps).*(1-du));
        F = F + (x <= eps+(M-eps).*du);
        F = F + (x <=-eps+(M+eps).*(1-dd));
        F = F + (x >= -eps+(m+eps).*dd);
        F = F + (m.*(du+dd) <= x <= M.*(dd+du));
        F = F + (0 <= z <=nx*mx);
        
        %F = F + (x > 0+(m-0).*(1-du));
        % F = F + (x < -0+(0+M).*(1-dd));
        %F = F + (0 <= z <=nx*mx);
        
      %  F = (m.*(du+dd) <= x <= M.*(du+dd)) + (sum([du dd],2) <= 1);
      %  F = F + (z == ((sum(du) + sum(dd))));
      %  F = F + (x > 1e-5+(-1+m).*(1-du));
      %  F = F + (x < -1e-5+(1+M).*(1-dd));
      %  F = F + (0 <= z <=nx*mx);
    otherwise
end


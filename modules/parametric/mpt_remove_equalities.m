function  Matrices = mpt_remove_equalities(Matrices,remove);
m  = size(Matrices.G,1);
nu = Matrices.nu;
nx = Matrices.nx;

Matrices.G(remove,:) = [];
Matrices.E(remove,:) = [];
Matrices.W(remove,:) = [];

% But add variable bounds (these where used to remove the rows)
Matrices.G = [Matrices.G;eye(nu);-eye(nu);zeros(2*nx,nu)];
Matrices.E = [Matrices.E;zeros(2*nu,nx);-eye(nx);eye(nx)];
Matrices.W = [Matrices.W;Matrices.ub(1:nu);-Matrices.lb(1:nu);Matrices.ub(nu+1:end);-Matrices.lb(nu+1:end)];
infbounds = find(isinf(Matrices.W) & (Matrices.W>0));
Matrices.G(infbounds,:) = [];
Matrices.E(infbounds,:) = [];
Matrices.W(infbounds,:) = [];

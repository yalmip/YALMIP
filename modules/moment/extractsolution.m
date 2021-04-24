function x_extract = extractsolution(momentdata,options)
%EXTRACTSOLUTIONS Tries to extract solutions from moment matrices
%
% 
%  xoptimal = extractsolution(momentstructure)
%
%   xoptimal    : Extracted solutions
%
%   momentdata : Problem data, obtained from SOLVEMOMENT (or SOLVESDP)
%   options    : Options structure from SDPSETTINGS
%
%   See also SOLVEMOMENT, SDPSETTINGS

moment = momentdata.moment;
x = momentdata.x;
monomials = momentdata.monomials;
n = momentdata.n;
d = momentdata.d;

[U,S,V,ranks] = numranks(moment);

if options.moment.extractrank>0
    % We try a extraction from highest order moment no matter what
    flat = length(moment);
    ranks(flat) = options.moment.extractrank;
else
    % Find a flat extension
    flat = d+min(find(ranks(1+d:end)-ranks(1:end-d)==0));
end

if ~isempty(flat)
    
    % Find a basis
    r = ranks(flat);
    V = U{flat}(:,1:r)*sqrt(diag(diag(S{flat}(1:r,1:r))));    
    if options.moment.rceftol >= 0
        [R,pivot] = rref(V',options.moment.rceftol);R = R';
    else
        % Try to find a reasonable tolerance by avoiding severly badly
        % conditioned R. Hack.., but seem to behave rather robustly.
        cV = cond(V);
        tol = 1e-10;
        [R,pivot] = rref(V',tol);R = R';
        while tol<1 & (cond(R)/cond(cV)>1e4)
            tol = tol*5;
            [R,pivot] = rref(V',tol);R = R';
        end
    end
    
    % Figure out multiplying matrices using YALMIP code
    w = monomials(pivot);
    for i = 1:n
        xw = x(i)*w;
        k = [];
        for j = 1:length(xw)
            k = [k;find(ismember(xw(j),monomials))];           
        end
        N{i} = R(k,:);
    end
    
    % Things missing in the basis...
    if ~all(cellfun('prodofsize',N)==length(w)^2)
       x_extract = {[]};
       return;
    end

    % Create random convex combination
    rands = rand(n,1);rands = rands/sum(rands);
    M = 0;
    for i = 1:n
        M = M + rands(i)*N{i};
    end

    [Q,T] = schur(M);
    % Extract solution
    for i = 1:r
        for j = 1:n
            x_extract{i}(j,1) =  Q(:,i)'*N{j}*Q(:,i);
        end
    end
    
    % Refine solutions v-Rw = e(x)=0
    if options.moment.refine>0
        xtemp = double(x);
        e = monomials(1:size(R,1))-R*w;        
        dedx = jacobian(e,x);
        for j = 1:r
            assign(x,x_extract{j});
            for i = 1:options.moment.refine
                assign(x,double(x)-double(dedx)\double(e));                
            end
            x_extract{j} = double(x);
        end
        assign(x,xtemp);
    end

else
    x_extract = {};
end

% Somewhat more stable rank detection
% looking for sharp drops in singular value
function [U,S,V,ranks] = numranks(moment)
for i = 1:length(moment)
    [U{i},S{i},V{i}] = svd(moment{i});
    s = diag(S{i});
    decay = s(2:end)./(eps+s(1:end-1));
    r = min(find(decay<1e-3));
    if isempty(r)
        ranks(i) = rank(moment{i},1e-8);
    else
        ranks(i) = min(r,rank(moment{i},1e-8));
    end
end


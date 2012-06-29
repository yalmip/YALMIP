function  test_pennon_vibration%(file, upper, rank, lambda_bar, resx, resy)


mbg_asserttrue(1);
return

draw = 1;

if nargin < 6
    draw = 0;
end

if nargin < 4
    lambda_bar = 0.0001;
end

if nargin < 2
    upper = 10;
end    
    
if nargin < 3
    rank = 0;
end

file = 'shape1';

load (file);


if upper > 0
    bounds = 1;
else
    bounds = 0;
end

tic;

compl = 5.0e-0; 
vol = nelem/3;

% the variables
u = sdpvar(nnod,1);
if rank > 0
    w = sdpvar(nnod,rank);
end

if bounds > 0
    r = sdpvar(nelem,1);
end

alpha = sdpvar(1,1);

lb = 1;
BIUALL = [];
NN = zeros(length(u));
for ie=1:nelem %for all elements
    len = Bdim(ie); %number of nonzeros in (Bi_1,...,Bi_nig)
    nb = len/nig; %number of nonzeros in Bi ...   
    SumBiuuBi = zeros(3,3);
    
    BIU = [];
    for ig=1:nig
        Bi = sparse(Brow(lb:lb+nb-1),Bcol(lb:lb+nb-1),Bval(lb:lb+nb-1),3,nnod);
        Biu = Bi*u;
        SumBiuuBi = SumBiuuBi + Biu*Biu';    
      
        if rank > 0
            Biw = Bi*w;
            SumBiuuBi = SumBiuuBi + Biw*Biw';    
            %for ir=1:rank
            %    ss = Bi*w(:,ir);
            %    SumBiuuBi = SumBiuuBi + ss*ss';
            %end
        end
        lb = lb + nb;
     %   BIU = [ BIU Biu];
    end
 
NN = NN | full(Bi)'*full(Bi);
    %SumBiuuBi = BIU*BIU';
    BuuB{ie} = 0.5*SumBiuuBi;
end

if rank > 0
    %Compue global mass matrix
    Mass = sparse(nnod,nnod);
    lb = 1;
    for ie = 1:nelem 
        len = Bdim(ie); %number of nonzeros in (Bi_1,...,Bi_nig)
        nb = len/nig; %number of nonzeros in Bi ...   
        elm=sparse(nnod,nnod);
        for ig=1:nig
            Bi = sparse(Brow(lb:lb+nb-1),Bcol(lb:lb+nb-1),Bval(lb:lb+nb-1),3,nnod);
            elm = elm + Bi'*Bi;  
            lb = lb + nb;
        end
        ss =  spones(diag(diag(elm)));
        MM{ie} = .5*ss;
    end
end

if rank > 0
    % This is to prevent rank-deficient solutions
    G = set(w(1,1)>=0);
    for ir=2:rank
    %    G = G+set(w(1,ir)<0);
        G = G+set(w(1,ir)>=0);
    end
    
    % Add matrix constraints
    for ie=1:nelem
        if bounds > 0
            G = G + set(BuuB{ie}-(lambda_bar*(trace(w'*MM{ie}*w))+alpha+r(ie))*eye(3)<=0);
        else
            G = G + set(BuuB{ie}-(lambda_bar*(trace(w'*MM{ie}*w))+alpha)*eye(3)<=0);
        end
    end

    %if rank > 1
    %    G = G + set(w'*w-0.01*eye(rank)>=0);
    %end
    
    if bounds > 0
        for ie=1:nelem
            G = G + set(r(ie)>=0);
        end
    end
else
    % Create 1-st constraint
    if bounds > 0
        G = set(BuuB{1}-(alpha+r(1))*eye(3)<=0);
        % Add constraint 2 to #elements
        for ie=2:nelem
            G = G + set(BuuB{ie}-(alpha+r(ie))*eye(3)<=0);
        end
        for ie=1:nelem
            G = G + set(r(ie)>=0);
        end
    else
        G = set(BuuB{1}-alpha*eye(3)<=0);
        % Add constraint 2 to #elements
        for ie=2:nelem
            G = G + set(BuuB{ie}-alpha*eye(3)<=0);
        end        
    end
    
end

% Create objective
if bounds > 0
    Obj = alpha*vol-F'*u+upper*r'*ones(nelem,1);
else
    Obj = alpha*vol-F'*u;
end

% Solve the problem
options = sdpsettings('solver','penbmi','usex0',1,'savedebug',0,...
    'penbmi.U0',1000,...
    'penbmi.P0',1000,...
    'penbmi.PEN_UP',.5,...
    'penbmi.ALPHA',1e-5,...
    'penbmi.P_EPS',1e-3,...%    'penbmi.PBM_EPS',1e-3,...
    'penbmi.PRECISION_2', 1e-4);  %1e-7

toc;
tic;

sol = solvesdp(G, Obj, options);
mbg_asserttrue(sol.problem == 0);




toc;

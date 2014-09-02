function [A,b] = generate_kernel_representation_data(N,N_unique,exponent_m2,exponent_p,p,options,p_base_parametric,ParametricIndicies,MonomIndicies,FirstRun)

persistent saveData

exponent_p_parametric = exponent_p(:,ParametricIndicies);
exponent_p_monoms = exponent_p(:,MonomIndicies);
pcoeffs = getbase(p);
if any(exponent_p_monoms(1,:))
    pcoeffs=pcoeffs(:,2:end); % No constant term in p
end
b = [];

parametric = full((~isempty(ParametricIndicies) & any(any(exponent_p_parametric))));

% For problems with a lot of similar cones, this saves some time
reuse = 0;
if ~isempty(saveData) && isequal(saveData.N,N) & ~FirstRun
    n = saveData.n;
    ind = saveData.ind;
    if  isequal(saveData.N_unique,N_unique) & isequal(saveData.exponent_m2,exponent_m2)% & isequal(saveData.epm,exponent_p_monoms)
        reuse = 1;
    end
else
    % Congruence partition sizes
    for k = 1:size(N,1)
        n(k) = size(N{k},1);
    end
    % Save old SOS definition
    saveData.N = N;
    saveData.n = n;
    saveData.N_unique = N_unique;
    saveData.exponent_m2 = exponent_m2;
    saveData.N_unique = N_unique;
end

if reuse & options.sos.reuse
    % Get old stuff
    if size(exponent_m2{1},2)==2 % Stupid (sos(parametric)) case
        ind = spalloc(1,1,0);
        ind(1)=1;
        allj = 1:size(exponent_p_monoms,1);
        used_in_p = ones(size(exponent_p_monoms,1),1);
    else
        allj = [];
        used_in_p = zeros(size(exponent_p_monoms,1),1);
        hash = randn(size(exponent_p_monoms,2),1);
        p_hash = exponent_p_monoms*hash;
        exponent_p_monoms_hash = exponent_p_monoms*hash;
        for i = 1:size(N_unique,1)
            monom = sparse(N_unique(i,3:end));
            j = find(exponent_p_monoms_hash == (monom*hash));
          
            if isempty(j)
                b = [b 0];
                allj(end+1,1) = 0;
            else
                used_in_p(j) = 1;
                allj(end+1,1:length(j)) = j(:)';
            end
        end
        ind = saveData.ind;
    end
else
    allj = [];
    used_in_p = zeros(size(exponent_p_monoms,1),1);
    if size(exponent_m2{1},2)==2 % Stupid (sos(parametric)) case
        ind = spalloc(1,1,0);
        ind(1)=1;
        allj = 1:size(exponent_p_monoms,1);
        used_in_p = ones(size(exponent_p_monoms,1),1);
    else
        % To speed up some searching, we random-hash data
        hash = randn(size(exponent_p_monoms,2),1);
        for k = 1:length(exponent_m2)
            if isempty(exponent_m2{k})
                exp_hash{k}=[];
            else
                exp_hash{k} = sparse((exponent_m2{k}(:,3:end)))*hash; % SPARSE NEEDED DUE TO STRANGE NUMERICS IN MATLAB ON 0s (the stuff will differ on last bit in hex format)
            end
        end

        p_hash = exponent_p_monoms*hash;
        ind = spalloc(size(N_unique,1),sum(n.^2),0);       

        for i = 1:size(N_unique,1)
            monom = N_unique(i,3:end);
            monom_hash = sparse(monom)*hash;
            LHS = 0;
            start = 0;
            for k = 1:size(N,1)
                j = find(exp_hash{k} == monom_hash);
                if ~isempty(j)
                    pos=exponent_m2{k}(j,1:2);
                    nss = pos(:,1);
                    mss = pos(:,2);
                    indicies = nss+(mss-1)*n(k);
                    ind(i,indicies+start) = ind(i,indicies+start) + 1;                                                                 
                end
                start = start + (n(k))^2;
                %                start = start + (matrixSOSsize*n(k))^2;
            end

            j = find(p_hash == monom_hash);

            if isempty(j)
                allj(end+1,1) = 0;
            else
                used_in_p(j) = 1;
                allj(end+1,1:length(j)) = j(:)';
            end
        end                     
    end
end
saveData.ind = ind;

% Some parametric terms in p(x,t) do not appear in v'Qv
% So these have to be added 0*Q = b
not_dealt_with  = find(used_in_p==0);
while ~isempty(not_dealt_with)
    %j = findrows(exponent_p_monoms,exponent_p_monoms(not_dealt_with(1),:));
    j = find(p_hash == p_hash(not_dealt_with(1)));
    allj(end+1,1:length(j)) = j(:)';
    used_in_p(j) = 1;
    not_dealt_with  = find(used_in_p==0);
    ind(end+1,1)=0;
end

matrixSOSsize = length(p);
if parametric
    % Inconsistent behaviour in MATLAB
    if size(allj,1)==1
        uu = [0;p_base_parametric];
        b = sum(uu(allj+1))';
    else
        b = [];
        for i = 1:matrixSOSsize
            for j = i:matrixSOSsize
                if i~=j
                    uu = [0;2*p_base_parametric(:,(i-1)*matrixSOSsize+j)];
                else
                    uu = [0;p_base_parametric(:,(i-1)*matrixSOSsize+j)];
                end
                b = [b sum(uu(allj+1),2)'];
            end
        end
    end
else
    if matrixSOSsize == 1
        uu = [zeros(size(pcoeffs,1),1) pcoeffs]';
        b = sum(uu(allj+1,:),2)';
    else
        b = [];
        for i = 1:matrixSOSsize
            for j = i:matrixSOSsize
                if i~=j
                    uu = [0;2*pcoeffs((i-1)*matrixSOSsize+j,:)'];
                else
                    uu = [0;pcoeffs((i-1)*matrixSOSsize+j,:)'];
                end
                b = [b;sum(uu(allj+1,:),2)'];
            end
        end
    end
    % uu = [0;pcoeffs(:)];
    % b = sum(uu(allj+1),2)';
end

b = b';
dualbase = ind;

j = 1;
A = cell(size(N,1),1);
for k = 1:size(N,1)
    if matrixSOSsize==1
        A{k} = dualbase(:,j:j+n(k)^2-1);
    else
        % Quick fix for matrix SOS case, should be optimized
        A{k} = inflate(dualbase(:,j:j+n(k)^2-1),matrixSOSsize,n(k));
    end
    j = j + n(k)^2;
end
b = b(:);



function newAi = inflate(Ai,matrixSOSsize,n);
% Quick fix for matrix SOS case, should be optimized
newAi = [];
newAi = [];
newAj = [];
newAk = [];
top = 1;
for i = 1:matrixSOSsize
    for r = i:matrixSOSsize
        for m = 1:size(Ai,1)
            ai = reshape(Ai(m,:),n,n);
            
            if 1
              %  V = spalloc(matrixSOSsize,matrixSOSsize,2);
              %  V(i,r)=1;
              %  V(r,i)=1;
              %  aii = kron(V,ai);
              %  aii = aii(:);
              %  [ii,jj,kk] = find(aii);
              % newAj = [newAj ii(:)'];
              %  newAi = [newAi repmat(top,1,length(ii))];
              %  newAk = [newAk kk(:)'];
                [dnewAj,dnewAi,dnewAk] = inflatelocal(ai,matrixSOSsize,r,i,top);
                newAj = [newAj dnewAj];
                newAi = [newAi dnewAi];
                newAk = [newAk dnewAk];
                % newAi = [newAi;ai(:)'];
            else
                 [dnewAjC,dnewAiC,dnewAkC] = inflatelocal(ai,matrixSOSsize,r,i,top);
                 [dnewAj,dnewAi,dnewAk] = inflatelocalnew(ai,matrixSOSsize,r,i,top,n);
                 AA=reshape(full(sparse(dnewAi*0+1,dnewAj,dnewAk,1,(matrixSOSsize*n)^2)),n*matrixSOSsize,[]);
                 AA2=reshape(full(sparse(dnewAiC*0+1,dnewAjC,dnewAkC,1,(matrixSOSsize*n)^2)),n*matrixSOSsize,[]);
                 
                 if norm(AA-AA2)>0
                     1
                 end
              
              
                newAj = [newAj dnewAj];
                newAi = [newAi dnewAi];
                newAk = [newAk dnewAk];
%                 
%                 [ii,jj,kk] = find(ai-diag(diag(ai)));
%                 iii = [(i-1)*n+ii;(r-1)*n+jj]
%                 jjj = [(r-1)*n+jj;(i-1)*n+ii]
%                 kkk = [kk;kk];
%                 indexi = repmat(top,1,length(iii));
%                 index = iii+(jjj-1)*matrixSOSsize*n
%                 newAj = [newAj index(:)'];
%                 newAi = [newAi indexi];
%                 newAk = [newAk kkk(:)'];
%                 
%                 [ii,jj,kk] = find(diag(diag(ai)));
%                 iii = [(i-1)*n+ii]
%                 jjj = [(r-1)*n+jj]
%                 kkk = [kk];
%                 indexi = repmat(top,1,length(iii));
%                 index = iii+(jjj-1)*matrixSOSsize*n
%                 newAj = [newAj index(:)'];
%                 newAi = [newAi indexi];
%                 newAk = [newAk kkk(:)'];
            end
            
            %sparse(indexi,index,kkk,1,(n*nA)^2)
            
            top = top+1;
        end
    end
end
newAi = sparse(newAi,newAj,newAk,top-1,(matrixSOSsize*n)^2);


function [Z,Q1,R] = sparsenull(A)

[Q,R] = qr(A');
n = max(find(sum(abs(R),2)));
Q1 = Q(:,1:n);
R = R(1:n,:);
Z = Q(:,n+1:end); % New basis


 function [dnewAj,dnewAi,dnewAk] = inflatelocal(ai,matrixSOSsize,r,i,top)
        V = spalloc(matrixSOSsize,matrixSOSsize,2);
        V(i,r)=1;
        V(r,i)=1;
        aii = kron(V,ai);
        aii = aii(:);
        [ii,jj,kk] = find(aii);
        
        dnewAj = ii(:)';
        dnewAi = repmat(top,1,length(ii));
        dnewAk = kk(:)';
        % newAi = [newAi;ai(:)'];
        
        
        
     function [dnewAj,dnewAi,dnewAk] = inflatelocalnew(ai,matrixSOSsize,r,i,top,n)
         
         if r==i
             [ii,jj,kk] = find(ai-diag(diag(ai)));
         else
             [ii,jj,kk] = find(ai);
         end               
         iii = [(i-1)*n+ii;(r-1)*n+jj];
         jjj = [(r-1)*n+jj;(i-1)*n+ii];
         kkk = [kk;kk];
         indexi = repmat(top,1,length(iii));
         index = iii+(jjj-1)*matrixSOSsize*n;
         dnewAj = [ index(:)'];
         dnewAi = [ indexi];
         dnewAk = [ kkk(:)'];
         
         [ii,jj,kk] = find(diag(diag(ai)));
         iii = [(i-1)*n+ii];
         jjj = [(r-1)*n+jj];
         kkk = [kk];
         indexi = repmat(top,1,length(iii));
         index = iii+(jjj-1)*matrixSOSsize*n;
         dnewAj = [dnewAj index(:)'];
         dnewAi = [dnewAi indexi];
         dnewAk = [dnewAk kkk(:)'];
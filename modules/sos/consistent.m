function keep = consistent(exponent_m,exponent_p)
% CONSISTENT Removes monomials using diagonal inconsistency
%
% V = CONSISTENT(V,P)
%
% Input
%  P : Scalar SDPVAR object
%  V : Vector with SDPVAR objects
%
% Output
%  V : Vector with SDPVAR objects
%
% Example:
%
% sdpvar x y
% p = 1+x^4*y^2+x^2*y^4;
% v = newtonmonoms(p);
% sdisplay(v)
% v = consistent(v,p);
% sdisplay(v)
%
% See also NEWTONREDUCE, NEWTONMONOMS, CONGRUENCEBLOCKS

% Author Johan Löfberg
% $Id: consistent.m,v 1.2 2009-11-24 14:53:54 joloef Exp $


% Default high level call with SDPVAR polynomials
% Convert to exponent form
sdpvarout = 0;
if isa(exponent_m,'sdpvar')
    z = depends(exponent_p);
    z = recover(unique([depends(exponent_p) depends(exponent_m)]));
    [exponent_p,p_base] = getexponentbase(exponent_p,z);
    exponent_m = getexponentbase(exponent_m,z);    
    sdpvarout = 1;
end

n = size(exponent_m,1);

% A bit silly, but I want to keep track 
% on the removed monoms, hence the output
% from this function is the indicies
% to the kept monoms
% 
% Append with indicies

indicies = (1:n)';

exponent_m = sparse(exponent_m);
exponent_p = sparse(exponent_p);
hash = rand(size(exponent_m,2),1);
exponent_m_hash = exponent_m*hash;
exponent_2m_hash = (2*exponent_m)*hash;
exponent_p_hash = (exponent_p)*hash;
% Check the actual candidates to be removed
candidates = [];
for i = 1:n
    m = 2*exponent_m(i,:);

%    index_in_p = findrows(exponent_p,m);
    index_in_p = findhash(exponent_p_hash,exponent_2m_hash(i),length(exponent_p_hash));
%    if ~isequal(index_in_p,index_in_p2)
%        [index_in_p index_in_p2]
%    end
    if isempty(index_in_p)
        candidates = [candidates;i];
    end
end

sums = full(sum(exponent_m,2));
for i = 0:max(sums)
    temp{i+1} = exponent_m(sums==i,:);  
end
maxsum = max(sums);

exponent_2m = 2*exponent_m;
exponent_m_transpose = exponent_m';
if ~isempty(candidates)
    removed = 1;    
    while removed 
        removelist = [];
        removed=0;
        n = size(exponent_m,1);
        possible = find(ismembc(indicies,candidates));
        j = 1;
        while j<=length(possible)
            i = possible(j);
            
            % Test this monom
           %m1 = exponent_m(i,:);
            
           %m_squared = 2*m1;
            m_squared = exponent_2m(i,:);
            % Can this be generated else-way, as mj+mk, j,k~=i
            terms = [];
            terms2 = [];
            m_combs = [1:1:i-1 i+1:1:n];
            m = length(m_combs);
            k = 1; 
                            
            k = 1; 
            %m_squared*hash-exponent_m(m_combs(k),:)*hash
            %full((m_squared*hash- - exponent_m(:,m_combs(k))'*hash)
            while (k<=m) & isempty(terms)
             %   m2 = exponent_m(m_combs(k),:);                
                m2 = exponent_m_transpose(:,m_combs(k))';                
                x = m_squared-m2;
%                x = 2*m1-m2;
                if all(x>=0)
                    terms = find(exponent_m_hash==x*hash);
                end               
                k = k+1;
            end
            
            if isempty(terms)
                removelist = [removelist i];
                removed = 1;
                j = inf;            
            end
            j = j+1;
        end
        exponent_m(removelist,:)=[];
        exponent_m_transpose(:,removelist)=[];        
        exponent_2m(removelist,:)=[];
        exponent_m_hash(removelist,:)=[];
        indicies(removelist)=[];        
    end
end
keep = indicies;

if sdpvarout
    keep = recovermonoms(exponent_m,z);
end

return

% ****************************
% MONOMIAL PRODUCTS
% ****************************
%  Faster (?) version
V = symminksum(exponent_m);


i = V(:,1)==V(:,2);
sumlist_diag = V(find(i),:);
sumlist = V;

% ******************************
removed = 1;
removelist = [];
slow = 1;
sumindex = sumlist(:,1:2);
sumdata = sumlist(:,3:end);

% % Partition data to make searches faster
merit = sum(sumdata,2) + sum(sumdata>0,2);
for i = min(merit):max(merit)
    part_sum_data{i+1}  = sumdata(find(merit==i),:);
    part_sum_index{i+1} = sumindex(find(merit==i),:);  
end

while removed
    removed = 0;
    removelist = [];
    for i = 1:size(sumlist_diag,1)
        this_monom = sumlist_diag(i,3:end);    
        index_in_p = findrows(exponent_p,this_monom);
        if isempty(index_in_p)
            m_this_monom = sum(this_monom) + sum(this_monom>0);            
            index_in_m = findrows(part_sum_data{m_this_monom+1},this_monom);            
            if ~isempty(index_in_m) & length(index_in_m)==1
                removed = 1;
                removelist = [removelist  sumlist_diag(i,1)];
            end
        end
    end
    
    for i = 1:length(part_sum_data)
        if ~isempty(part_sum_data{i})
            index = find(~(ismember(part_sum_index{i}(:,1),removelist) | ismember(part_sum_index{i}(:,2),removelist)));
            part_sum_index{i} = part_sum_index{i}(index,:);
            part_sum_data{i}  = part_sum_data{i}(index,:);
        end
    end 
end
keep = [];
for i = 1:length(part_sum_data)
    if ~isempty(part_sum_data{i})
        keep = [keep;part_sum_index{i}];
    end
end 
keep = unique(keep);
%keep = unique(sumindex);

if sdpvarout
    keep = recovermonoms(exponent_m,z);
end






function V = symminksum(exponent_m)

n = size(exponent_m,1);

%H1 = kron(speye(n),ones(n,1));
%H2 = kron(ones(n,1),speye(n));
%V0 = H1*exponent_m + repmat(exponent_m,n,1);
V0 = kron(exponent_m,ones(n,1)) + repmat(exponent_m,n,1);

ind = [];
k = 0;
indicies = [];
for i = 1:size(exponent_m,1)
    ind = [ind;k+(i:n)'];k = k+n;
    indicies =[indicies;(i:n)' repmat(i,n-i+1,1)];
end
V0 = V0(ind,:);
V = [indicies V0];

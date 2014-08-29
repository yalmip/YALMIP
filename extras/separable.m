function  exponent_m = separable(exponent_m,exponent_p,options);
%SEPARABLE Internal function, not used

% %exponent_m(sum((exponent_m>0),2)>2,:)=[];
% 
% card = max(sum((exponent_p>0),2));
% 
% n_less = exponent_m(sum((exponent_m>0),2)<card,:);
% n_equ  = exponent_m(sum((exponent_m>0),2)==card,:);
% n_larg = exponent_m(sum((exponent_m>0),2)>card,:);
% 
% A = minksum(n_less,n_less);
% B = minksum(n_less,n_equ);
% C = minksum(n_less,n_larg);
% D = minksum(n_equ,n_equ);
% E = minksum(n_equ,n_larg);
% F = minksum(n_larg,n_larg);

disconnected = [];
for i = 1:size(exponent_p,2)
    for j = i+1:size(exponent_p,2)
        if ~any(exponent_p(:,i) & exponent_p(:,j))
            disconnected = [disconnected;i j];
        end
    end
end

for i = 1:size(disconnected,1)
    j = disconnected(i,1);
    k = disconnected(i,2);
    n0 = find(~exponent_m(:,j) & ~exponent_m(:,k));
    nx = find(exponent_m(:,j) & ~exponent_m(:,k));
    nz = find(~exponent_m(:,j) & exponent_m(:,k));
    nxz = find(exponent_m(:,j) & exponent_m(:,k));
    
%     m0 = exponent_m(n0,:);
%     mx = exponent_m(nx,:);
%     mz = exponent_m(nz,:);
%     mxz = exponent_m(nxz,:);
%     
%     from_E = minksum(mx,mz);
%     from_B = minksum([m0;mx;mz],mxz);
%     from_C = minksum(mxz,mxz); 
%     m_e = exponent_m(union(nx,nz),:)
%     m_cb = exponent_m(union(nx,nz),:)
    exponent_m = exponent_m([n0;nx;nz],:);
end


function msum = minksum(a,b);
msum = [];
for i = 1:size(a,1)
    for j = i:size(b,1)
        msum = [msum;a(i,:)+b(j,:)];
    end
end
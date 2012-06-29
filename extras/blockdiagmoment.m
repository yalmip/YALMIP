function M = blockdiagmoment(obj,M)

rt = [];

r_old = [];
go_on = 1;
while go_on
    exponent_p = getexponentbase(obj+sum(sum(rt)),recover(depends(obj)));
    Msparsity = zeros(size(M{end},1));
    vars = recover(depends(obj));
    for i = 1:size(M{end},1)
        for j = i:size(M{end},1)
            exponent_Mij = getexponentbase(M{end}(i,j),vars);
            if ~isempty(findrows(exponent_p,exponent_Mij)) & sum(exponent_Mij)~=0
                Msparsity(i,j) = 1;
                Msparsity(j,i) = 1;
            end
        end
    end

    [p,q,r,s] = dmperm(Msparsity+eye(length(Msparsity)));
    if isequal(r,r_old)
        go_on = 0;
    else
        Mblocked = [];
        Mrs = M{end}(p,p);       
        MM = ones(length(M{end}));
        blocks = zeros(1,length(M{end}));
        Ms = Msparsity(p,p);
        for i = 1:length(r)-1
            blocks(r(i+1)-r(i)) = blocks(r(i+1)-r(i)) + 1;            
            Mtest = Ms(r(i):(r(i+1)-1),r(i):(r(i+1)-1));
            if nnz(Mtest)==0
            else                
                MM = blkdiag(MM,ones(r(i+1)-r(i)));
                Mblocked = blkdiag(Mblocked,Mrs(r(i):(r(i+1)-1),r(i):(r(i+1)-1)));
            end
        end
        string = 'Blocks : ';
        for i = 1:length(blocks)
            if blocks(i)>0
                string = [string num2str(i) 'x' num2str(i) '(' num2str(blocks(i)) ') ' ];
            end
        end
        disp(string)
        rt = Mblocked(:);
        r_old = r;
    end
end
M{end} = Mblocked;

% 
% mt = yalmip('monomtable');
% hash = randn(size(mt,2),1);
% for i = 1:length(u{end})
%     for j = 1:length(u{end})
%         if i==1
%             v1 = zeros(1,size(mt,2));
%         else
%             vi = getvariables(u{end}(i));
%             v1 = mt(vi,:);
%         end
%         if j==1
%             v2 = zeros(1,size(mt,2));
%         else
%             vj = getvariables(u{end}(j));
%             v2 = mt(vj,:);
%         end
%         v = v1+v2;
%         vh = v*hash;
%         uiuj(i,j) = vh;
%     end
% end


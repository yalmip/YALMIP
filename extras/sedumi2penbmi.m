function pen = sedumi2penbmi(F_struc,c,Q,K,monomtable,options,x0)

% Linear before bilinearize
temp = sum(monomtable,2)>1;
tempnonlinearindicies = find(temp);
templinearindicies = find(~temp);

if isempty(monomtable)
    monomtable = eye(length(c));
end

temp = sum(monomtable,2)>1;
nonlinearindicies = find(temp);
linearindicies = find(~temp);
c0 = c;
c  = c(linearindicies);
Q = Q(linearindicies,linearindicies);
nonlinear_scalars = [];

% Any non-linear scalar inequalities?
% Move these to the BMI part
if K.l>0
    nonlinear_scalars = find(any(full(F_struc(1:K.l,[nonlinearindicies(:)'+1])),2));
    if ~isempty(nonlinear_scalars)
        Kold = K;
        linear_scalars = setdiff(1:K.l,nonlinear_scalars);
        F_struc = [F_struc(linear_scalars,:);F_struc(nonlinear_scalars,:);F_struc(K.l+1:end,:)];
        K.l = K.l-length(nonlinear_scalars);
        if (length(K.s)==1) & (K.s==0)
            K.s    = [repmat(1,1,length(nonlinear_scalars))];
            K.rank = repmat(1,1,length(nonlinear_scalars)); 
        else
            K.s    = [repmat(1,1,length(nonlinear_scalars)) K.s];
            K.rank = [repmat(1,1,length(nonlinear_scalars)) K.rank];
        end
    end
end

if ~isempty(F_struc)
    pen = sedumi2pen(F_struc(:,[1 linearindicies(:)'+1]),K,c,x0);
else
    pen = sedumi2pen([],K,c,x0);
end

if ~isempty(nonlinearindicies)
    bmi = sedumi2pen(F_struc(:,[nonlinearindicies(:)'+1]),K,[],[]);
    pen.ki_dim = bmi.ai_dim;
    % Nonlinear index
    pen.ki_dim = bmi.ai_dim;
    pen.ki_row = bmi.ai_row;
    pen.ki_col = bmi.ai_col;
    pen.ki_nzs = bmi.ai_nzs;
    pen.ki_val = bmi.ai_val;
    for i = 1:length(bmi.ai_idx)
        nl = nonlinearindicies(1+bmi.ai_idx(i));
        v = find(monomtable(nl,:));
        if length(v)==1
            v(2)=v(1);
        end      
        pen.ki_idx(i)=find(linearindicies == v(1));
        pen.kj_idx(i)=find(linearindicies == v(2));
    end
else
    pen.ki_dim = 0*pen.ai_dim;
    pen.ki_row = 0;
    pen.ki_col = 0;
    pen.ki_nzs = 0;
    pen.ki_idx = 0;
    pen.kj_idx = 0;
    pen.kj_val = 0;    
end

if nnz(Q)>0
    [row,col,vals] = find(triu(Q));
    pen.q_nzs = length(row);
    pen.Q_nzs = length(row);%Bug in PENLABs interface?
    pen.q_val = vals';
    pen.q_col = col'-1;
    pen.q_row = row'-1;
else
    pen.q_nzs = 0;
    pen.q_val = 0;
    pen.q_col = 0;
    pen.q_row = 0;    
end

ops = struct2cell(options.penbmi);ops = [ops{1:end}];
pen.ioptions = ops(1:12);
pen.foptions = ops(13:end);
pen.ioptions(4) = max(0,min(3,options.verbose+1));
if pen.ioptions(4)==1
    pen.ioptions(4)=0;
end

% ****************************************
% UNCOMMENT THIS FOR PENBMI version 1
% ****************************************
%pen.ioptions = pen.ioptions(1:8);
%pen.foptions = pen.foptions(1:8);

if ~isempty(x0)    
    pen.x0(isnan(pen.x0)) = 0;
    pen.x0 = x0(linearindicies);    
    pen.x0 = pen.x0(:)';
end
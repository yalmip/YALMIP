function pen = sedumi2pen(F_struc,K,c,x0)
%SEDUMI2MPEN Internal function to convert SeDuMi structure to format needed in PENNON

% General data
pen.vars = length(c);
pen.constr = K.l;
if length(K.s)>1
    pen.mconstr = length(K.s);
else
    if K.s==0
        pen.mconstr=0;
    else
        pen.mconstr=1;
    end
end
pen.msizes = K.s;
pen.fobj = c(:)';
pen.x0 = zeros(1,length(c));

% Linear constraints
if K.l>0
    pen.ci = full(F_struc(K.f+1:K.f+K.l,1))';
    pen.bi_dim = zeros(1,K.f+K.l);
    [variables,constr,val] = find(-F_struc(1:K.l,2:end)');
    pen.bi_idx = variables(:)'-1;
    pen.bi_val = val(:)';
    j = 1;

    if ~isempty(constr)
        rr = histc(constr,[1:max(constr)]);
        uu = find(rr);
        pen.bi_dim(uu) = rr(uu);
        %uns = uniquestripped(constr);
        %for j = 1:length(uns)
        %    pen.bi_dim(uns(j)) = sum(constr==uns(j));
        %end
    end

else
    pen.ci = 0;
    pen.bi_dim = 0;
    pen.bi_idx = 0;
    pen.bi_val = 0;
end

% Semidefnite constraints
if K.s>0
    top = K.l+K.f+1;
    constraints = 1;
    pen.ai_dim = zeros(1,length(K.s));

    % First, optimized code for scalar SDPs (happens for nonlinear scalar constraints)
    if all(K.s==1)

        Avec = -F_struc(top:top+length(K.s)-1,:);
        [variables,constr,val] = find(Avec');
        j = 1;
        n_nz = length(val);
        pen.ai_nzs = repmat(1,1,n_nz);
        pen.ai_row = repmat(0,1,n_nz);
        pen.ai_col = pen.ai_row;%repmat(0,1,n_nz);
        pen.ai_val = val(:)';
        pen.ai_idx = variables(:)'-1;
        while j<=length(val)
            first_constraint = constr(j);
            n_nz = 1;
            while j+n_nz<=length(val) & constr(j+n_nz)==first_constraint
                n_nz = n_nz + 1;
            end
            pen.ai_dim(first_constraint) = n_nz;
            j = j + n_nz;
        end
    else
        pen.ai_idx = [];
        pen.ai_val = [];
        pen.ai_row = [];
        pen.ai_col = [];
        pen.ai_nzs = [];        
        while (constraints<=length(K.s))
            n = K.s(constraints);           
            Avec = -F_struc(top:top+n^2-1,:);
            picker = reshape(1:n^2,n,n);picker = tril(picker-diag(diag(picker)));picker = find(picker(:));            
            [rowcols,varindicies,vals]=find(Avec);
            %the_col = 1+floor((rowcols-1)/n);
            %the_row = rowcols-(the_col-1)*n;
            %removethese = the_row > the_col;
            removethese = find(ismembc(rowcols,picker));
            rowcols(removethese) =[];
            varindicies(removethese)=[];
            vals(removethese)=[];
            if ~isempty(rowcols)
                cols = ceil(rowcols/n);
                rows = rowcols - n*(cols-1);
                uns = uniquestripped(varindicies);
               
                difference = diff([varindicies(:) ; max(varindicies)+1]);
                count = diff(find([1 ; difference]))';
                pen.ai_nzs = [pen.ai_nzs count];
                                
                pen.ai_dim(constraints) =  length(uniquestripped(varindicies));
                pen.ai_idx = [pen.ai_idx uns(:)'-1];
                pen.ai_row = [pen.ai_row rows(:)'-1];
                pen.ai_col = [pen.ai_col cols(:)'-1];
                pen.ai_val = [pen.ai_val vals(:)'];
            end
            constraints = constraints+1;
            top = top+n*n;
        end
    end
else
    pen.ai_dim = 0;
    pen.ai_idx = 0;
    pen.ai_val = 0;
    pen.ai_row = 0;
    pen.ai_col = 0;
    pen.ai_nzs = 0;
end
function pcut = addBilinearVariableCuts(p)

pcut = p;

z = p.bilinears(:,1);
x = p.bilinears(:,2);
y = p.bilinears(:,3);

nn = length(p.c);

still_uncertain = find(abs(p.lb(z)-p.ub(z))>1e-8);
if ~isempty(still_uncertain)

    x_lb = p.lb(x);
    x_ub = p.ub(x);
    y_lb = p.lb(y);
    y_ub = p.ub(y);
    m = length(x);
    one = ones(m,1);
    general_vals =[x_lb.*y_lb one -y_lb -x_lb,x_ub.*y_ub one -y_ub -x_ub,-x_ub.*y_lb -one y_lb x_ub,-x_lb.*y_ub -one y_ub x_lb]';
    general_cols = [one z+1 x+1 y+1 one z+1  x+1 y+1 one z+1 x+1 y+1 one z+1 x+1 y+1]';
    general_row = [1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];

    if 0
        quadratic_row = [1;1;1;2;2 ;2; 3; 3; 3];
        quadratic_cols =  [one  z+1 x+1 one z+1 x+1 one  z+1 x+1]';
        quadratic_vals = [-x_ub.*x_lb -one x_lb+x_ub x_lb.*y_lb one -y_lb-x_lb x_ub.*y_ub one -y_ub-x_ub]';
    else
        quadratic_row = [1;1;1;2;2 ;2;3; 3; 3;4;4;4;5;5;5;6;6;6];
        quadratic_cols =  [one  z+1 x+1 ,one z+1 x+1 ,one  z+1 x+1, one z+1 x+1 ,one z+1 x+1, one z+1 x+1]';
        x1 = (3*x_ub+x_lb)/4;
        x2 = (x_ub+3*x_lb)/4;
        x3 = (x_ub+x_lb)/2;
        quadratic_vals = [-x_ub.*x_lb -one x_lb+x_ub x_lb.*y_lb one -y_lb-x_lb x_ub.*y_ub one -y_ub-x_ub  x1.*x1 one -x1-x1 x2.*x2 one -x2-x2 x3.*x3 one -x3-x3]';
    end
    m = 1+length(p.c);
    rows = [];
    cols = [];
    vals = [];
    nrow = 0;

    for i =still_uncertain(:)'
        x = p.bilinears(i,2);
        y = p.bilinears(i,3);
        if x~=y
            rows = [rows;general_row+nrow];
            vals = [vals;general_vals(:,i)];
            cols = [cols;general_cols(:,i)];
            nrow = nrow + 4;
        else
            col = quadratic_cols(:,i);
            val = quadratic_vals(:,i);

            rows = [rows;quadratic_row+nrow];
            vals = [vals;val];
            cols = [cols;col];
            nrow = nrow + max(quadratic_row);
        end
    end

    F_temp = sparse(rows,cols,vals,nrow,m);   
    keep = find(~isinf(F_temp(:,1)) & ~isnan(F_temp(:,1)));
    F_temp = F_temp(keep,:);
    pcut.F_struc = [F_temp;pcut.F_struc];
    pcut.K.l = pcut.K.l+size(F_temp,1);
end


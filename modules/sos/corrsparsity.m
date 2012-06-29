function [C,D] = corrsparsity(exponent_p_monoms,options);

if options.sos.csp
    n = size(exponent_p_monoms,2);
    C = zeros(n,n);
    for i = 1:size(exponent_p_monoms,1)
        [col,row] = find(exponent_p_monoms(i,:));
        for j = 1:length(row)
            C(row(j),row(j)) = 1;
            for k = 2:length(row)
                C(row(j),row(k)) = 1;
                C(row(k),row(j)) = 1;
            end
        end
    end

    R = chol(C+100*eye(n));
    for i = 1:size(R,1)
        [row,col] = find(R(i,:));
        if i>1
            is_in = 0;
            for j = 1:length(D)
                if all(ismember(col,D{j}))
                    is_in = 1;
                end
            end
            if ~is_in
                D{end+1} = col;
            end
        else
            D{1} = col;
        end
    end
    if length(D)>1 &  options.verbose>0
        the_text = 'Detecting correlative sparsity..';

        [uu,ii,oo] = unique(cellfun('prodofsize',D));
        for i = 1:length(uu)
            n_this = length(find(oo==i));

            the_text = [the_text num2str(uu(i)) '(' num2str(n_this) ')' ' '];
        end       
        disp(the_text);
    end
else
    C = ones(size(exponent_p_monoms,2));
    D{1} = 1:size(exponent_p_monoms,2);
end

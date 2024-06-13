function N = blockmonomialreduction(exponent_p,N,options)
%BLOCKMONOMIALREDUCTION Internal function to reduce monomials in SOS problem

if options.sos.inconsistent & ~options.sos.csp
    n_removed = 1;
    while n_removed>0
        t = cputime;
        N_unique = [];
        for i = 1:length(N)
            for j = 1:size(N{i},1)
                for k = j:size(N{i},1)
                    N_unique = [N_unique;N{i}(j,:)+N{i}(k,:)];
                end
            end
        end
        
        n_removed = 0;
        for i = 1:length(N)
            rmv = [];
            for j = 1:size(N{i},1)
                nn = N{i}(j,:)*2;
                if isempty(findrows(exponent_p,nn)) & (length(findrows(N_unique,nn))==1)
                    rmv = [rmv j];                
                end
            end
            N{i}(rmv,:)=[];
            n_removed = n_removed + length(rmv);
        end
        M=[];k=1;
        for i = 1:length(N)
            if ~isempty(N{i})
                M{k,1}=N{i};k=k+1;
            end
        end
        N=M;
        t = cputime-t;
        if n_removed>0
            the_text = 'Block-diagonal inconsistensies..Keeping ';
            for i = 1:length(N)
                the_text = [the_text num2str(size(N{i},1)) 'x' num2str(size(N{i},1)) ' '];
            end
            the_text = [the_text ' (' num2str(t) 'sec)']; 
            disp(the_text);
        end
    end
end

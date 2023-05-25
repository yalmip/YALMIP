function cliques = derivecliques(a,b)

a = a(:)';
[a_,index] = sort(a,'ascend');
if a_(end-1)+a_(end) <= b
    % There cannot be any cliques here
    cliques = {};
else
    j = 1;
    n = length(a);
    % Search for first violation
    % FIXME should be done with binary search
    while 1
        if a_(j)+a_(j+1) > b
            break
        end
        j = j+1;
    end
    cliquestart = j;
    cliques{1} = sort(index(cliquestart:end));    
    % Brute-froce code to derive cliques involving the rest of the
    % variables
    for outside = 1:j-1
        for k = cliquestart:n
            if a_(outside) + a_(k) > b
                cliques{end+1} = sort([index(outside) index(k:end)]);
                break
            end
        end
    end
end

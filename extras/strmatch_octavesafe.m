function idx = strmatch_octavesafe(S, A, exact)
    if nargin < 3
        idx = find(strncmp(S, A, numel(S)));
    else
        % Exact match
        idx = find(strcmp(S, A));
    end
end


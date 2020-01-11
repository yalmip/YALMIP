function idx = strmatch_octavesafe(S, A)
    idx = find(strncmp(S, A, numel(S)));
end


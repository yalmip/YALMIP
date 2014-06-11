function members=ismembcYALMIP(a,b)

% ismembc is fast, but does not exist in octave
try
    members = ismembc(a,b);
catch
    members = ismember(a,b);
end

  
  

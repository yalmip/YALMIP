function indicies = find(x)
base = x.basis;
vars = x.lmi_variables;
indicies = find(any(base,2));

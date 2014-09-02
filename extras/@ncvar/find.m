function indicies = find(x)
%FIND (overloaded)

base = x.basis;
vars = x.lmi_variables;
indicies = find(any(base,2));

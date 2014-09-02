function Z = conj(Y)
%CONJ (overloaded)

Z = Y;
Z.basis = conj(Z.basis);
% Reset info about conic terms
Z.conicinfo = [0 0];
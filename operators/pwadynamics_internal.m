function [F,vars] = pwadynamics_internal(first_coord,varargin)

% Hack to figure out all the sorted variables, not just this index.
% Sort is implemented in a slighlt different way (general feature in
% future versions) that allows one element in an operator to modell all
% elements. Reduces the number of calls to the operator code.

n = length(varargin{1});
xplus = recover(getvariables(first_coord):getvariables(first_coord)+n-1);
vars = getvariables(xplus);
F = set([]);
for i = 1:(length(varargin)/2)
    fi{i} = varargin{2*i-1};
    R{i} =  varargin{2*i} + set(xplus == fi{i});
    F = F + R{i};
end

[F,t] = hull(R{:});F = F + set(binary(t));


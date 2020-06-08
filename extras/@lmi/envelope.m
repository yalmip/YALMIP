function [E,P] = envelope(C,x)
%ENVELOPE Create linear approximation of envelope of nonlinear constraint
%
%     [E,P] = envelope(C,x)
%
% C: Constraint object involving nonlinear expresions
% x: Optional: The linear variables of interest (projection onto these will be performed)
% E: Constraint object representing the envelope approximation
% P: Optional: Polyhedral object representing E (MPT Toolbox object)
%              (requires a second argument with variables to project
%              envelope to)
%
% Examples
%
% In order to derive the outer approximation in the original linear
% variables, the model has to be projected onto the linear variables
% This projetion based format requires MPT3
%   sdpvar x u
%   E = envelope([-1 <= x <= 1, u == x^3],[x;u]);
%   plot(E)
%   xx = (-1:0.01:1);hold on;plot(xx,xx.^3)
%
%   E = envelope([-1 <= x <= 1,x+sin(pi*x) <= u <= 4-x^2],[x;u]);
%   plot(E)
%   xx = (-1:0.01:1);hold on;plot(xx,xx+sin(pi*xx),xx,4-xx.^2)
%
% A polyhedral object in MPT format representing E can be obtained as a
% second output.
%
%   [E,P] = envelope([-1 <= x <= 1,x+sin(pi*x) <= u <= 4-x^2],[x;u]);
%
% Alternatively, to avoid an expensive projection,  we can create a model
% which adds the outer approximation cuts for the envelopes of the
% nonlinear variables, but keep the nonlinear variables, and then plot the
% projection of the outer approximation, keeping in mind that we now have
% to relax the nonlinear variables (this is the model which thus would be
% used in a branch&bound scheme) 
%   E = envelope([-1 <= x <= 1,x+sin(pi*x) <= u <= 4-x^2]);
%   plot(E,[x;u],[],[],sdpsettings('relax',1))
%   xx = (-1:0.01:1);hold on;plot(xx,xx+sin(pi*xx),xx,4-xx.^2)

% Author Johan Löfberg

[aux1,aux2,aux3,p] = export(C,[],sdpsettings('solver','bmibnb'));

if isempty(p)
    if ~isempty(aux3)
        aux3.info
    end
    error('Failed to export a model')
end

% Copied from bmibnb
p.EqualityConstraintState = ones(p.K.f,1);
p.InequalityConstraintState = ones(p.K.l,1);
p = compile_nonlinear_table(p);
%p = propagatequadratics(p,inf,-inf);
p = propagate_bounds_from_arbitrary_quadratics(p);
p.high_monom_model=[];
p.originalModel = p;
p = presolveOneMagicRound(p);   
[p,changed] = convert_sigmonial_to_sdpfun(p);
if changed
    p = compile_nonlinear_table(p);
end
[p,changed] = convert_polynomial_to_quadratic(p);
if changed
    p = compile_nonlinear_table(p);
    p.EqualityConstraintState = ones(p.K.f,1);
    p.InequalityConstraintState = ones(p.K.l,1);
end
p = presolveOneMagicRound(p);  

% Copied from solvelower
p_cut = p;
p_cut = addNormBoundCut(p_cut);
p_cut = addBilinearVariableCuts(p_cut);
p_cut = addEvalVariableCuts(p_cut);
p_cut = addMonomialCuts(p_cut);
p_cut = addMonomialTowerCuts(p_cut);
p_cut = addSinCosCuts(p_cut);

p_cut = mergeBoundsToModel(p_cut);
if nargin > 1
    % Now project onto the variables of interest
    for i = 1:length(x)
        xi(i) = find(getvariables(x(i)) == p.used_variables);
    end
    A = -p_cut.F_struc(:,2:end);
    b = p_cut.F_struc(:,1);
    
    Akeep = A(:,xi);
    A(:,xi)=[];
    
    A = [Akeep A];
    
    Ae = A(1:p_cut.K.f,:);
    be = b(1:p_cut.K.f,:);
    A = A(1+p_cut.K.f:end,:);
    b = b(1+p_cut.K.f:end,:);
    P = Polyhedron('A',A,'b',b,'Ae',full(Ae),'be',full(be));
    P = projection(P,1:length(xi));
    E = ismember(x,P);
else
    z = recover(p_cut.used_variables);
    % We might have introduced some local modelling variables here
    m = size(p.F_struc,2)-1-length(z);
    if m>0
        z = [z;sdpvar(m,1)];
    end
    E = [];
    top = 1;
    if p_cut.K.f > 0
        E = [E,p_cut.F_struc(1:p_cut.K.f,:)*[1;z]==0];
    end
    if p_cut.K.l > 0
        E = [E,p_cut.F_struc(1+p_cut.K.f:p_cut.K.f + p_cut.K.l,:)*[1;z]>=0];
    end
    if p_cut.K.q > 0
        top = 1 + p_cut.K.f + p_cut.K.l ;
        for i = 1:length(p_cut.K.q)
            n = p_cut.K.q(i);
            M = p_cut.F_struc(top:top+n-1,:)*[1;z];
            E = [E, cone(M)];
            top = top + n;
        end
    end
    if p_cut.K.s(1) > 0
        top = 1 + p_cut.K.f + p_cut.K.l + sum(p_cut.K.q);
        for i = 1:length(p_cut.K.s)
            n = p_cut.K.s(i);
            M = p_cut.F_struc(top:top+n^2-1,:)*[1;z];
            top = top + n^2;
            E = [E,reshape(M,n,n)>=0];
        end
    end 
end

function p = mergeBoundsToModel(p);

A = [];
b = [];
if ~isempty(p.lb)
    A = [eye(length(p.c))];
    b = p.ub;
end
if ~isempty(p.ub)
    A = [A;-eye(length(p.c))];
    b = [b;-p.lb];
end
infbounds = find(isinf(b));
A(infbounds,:)=[];
b(infbounds)=[];
if length(b)>0
    p.F_struc = [p.F_struc(1:p.K.f,:);[b -A];p.F_struc(p.K.f+1:end,:)];
    p.K.l = p.K.l + length(b);
end








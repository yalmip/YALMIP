function [sol,x_extract,momentsstructure,sosout,Fnew,obj] = solvemoment(F,obj,options,k)
%SOLVEMOMENT Application of Lasserre's moment-method for polynomial programming
%
%   min        h(x)
%   subject to
%              F(x) >= 0,
%
%   [DIAGNOSTIC,X,MOMENT,SOS,Flinear,Objlinear] = SOLVEMOMENT(F,h,options,k)
%
%   diagnostic : Struct with diagnostics
%   x          : Extracted global solutions
%   moment     : Structure with moments, various variables needed to recover solution
%   sos        : SOS decomposition {max t s.t h-t = p0+sum(pi*Fi), pi = vi'*Qi*vi}
%   Flinear    : The linearized constraints
%   Objlinear  : The linearized objective
%
%   Input
%      F       : SET object with polynomial inequalities and equalities.
%      h       : SDPVAR object describing the polynomial h(x).
%      options : solver options from SDPSETTINGS.
%      k       : Level of relaxation. If empty or not given, smallest possible applied.
%
%   The behaviour of the moment relaxation can be controlled
%   using the fields 'moment' in SDPSETTINGS
%
%      moment.refine      : Perform #refine Newton iterations in extracation of global solutions.
%                           This can improve numerical accuracy of extracted solutions in some cases.
%      moment.extractrank : Try (forcefully) to extract #extractrank global solutions.
%                           This feature should normally not be used and is default 0.
%      moment.rceftol     : Tolerance during Gaussian elimination used in extraction of global solutions.
%                           Default is -1 which means heuristic choice by YALMIP.
%
%   Some of the fields are only used when the moment relaxation is called
%   indirectly from SOLVESDP.
%
%      moment.solver : SDP solver used in moment relxation. Default ''
%      moment.order  : Order of relxation. Default [] meaning lowest possible.
%
%   See also SDPVAR, SET, SDPSETTINGS, SOLVESDP

% Author Johan Löfberg, Philipp Rostalski  Update
% $Id: solvemoment.m,v 1.8 2007/05/28 09:09:42 joloef Exp $
%
% Updated equality constraint handling, 2010/08/02

if nargin ==0
    help solvemoment
    return
end

if nargin<2
    obj=[];
end

if (nargin>=3) & (isa(options,'double') & ~isempty(options))
    help solvemoment
    error('Order of arguments have changed in solvemoment. Update code');
end

if nargin<3 | (isempty(options))
    options = sdpsettings;
end

if strcmp(options.solver,'sparsepop')
    if nargin >= 4        
        sol = solvesdp(F,obj,sdpsettings(options,'sparsepop.relaxorder',k));
    else
        sol = solvesdp(F,obj,options);
    end
    return
end

% Relaxation-order given?
if nargin<4
    k = options.moment.order;
end

% Check for wrong syntax
if ~isempty(F) & ~isa(F,'lmi') & ~isa(F,'constraint')
    error('First argument should be a SET object')
end

if isa(F,'constraint')
    F = lmi(F);
end

% Take care of rational expressions
[F,failure] = expandmodel(F,obj);
if failure
    error('Could not expand model (rational functions, min/max etc). Avoid nonlinear operators in moment problems.');
end

[Fnew,obj,M,k,x,u,n,deg,linears,nonlinears,vecConstraints,isinequality,ulong] = momentmodel(F,obj,k,1);

% No objective, minimize trace on moment-matrix instead
if isempty(obj)
    obj = trace(M{k+1});
end

% Solve
sol = solvesdp(Fnew,obj,sdpsettings(options,'relax',1));
assign(nonlinears,value(linears));

% Construct SOS decompositions if the user wants these
if nargout >= 4
    sosout.t = relaxdouble(obj);
    sosout.Q0 = dual(Fnew(1));
    sosout.v0 = u{end};
    sosout.p0 = u{end}'*dual(Fnew(1))*u{end};
    for i = 1:length(vecConstraints)
        if isinequality(i)
            sosout.Qi{i} = dual(Fnew(i+1));
            sosout.vi{i} = u{end}(1:length(sosout.Qi{i}));
            sosout.pi{i} = sosout.vi{i}'*sosout.Qi{i}*sosout.vi{i};
        else
            sosout.Qi{i} = dual(Fnew(i+1));
            sosout.vi{i} = ulong{end}(1:length(sosout.Qi{i}));
            sosout.pi{i} = sosout.Qi{i}'*sosout.vi{i};
        end
    end
end

% Get the moment matrices
% M{end} = M_original;
for i = 1:k+1
    moments{i} = relaxdouble(M{i});
end

% Extract solutions if possible (at-least fesible and unbounded)
momentsstructure.moment  = moments;
momentsstructure.x       = x;
momentsstructure.monomials = u{k};
momentsstructure.n       = n;
momentsstructure.d       = max(1,ceil(max(deg)/2));
x_extract = {};
if nargout>=2 & ~(sol.problem == 1 | sol.problem == 2)
    momentsstructure.moment  = moments;
    momentsstructure.x       = x;
    momentsstructure.monomials = u{k};
    momentsstructure.n       = n;
    momentsstructure.d       = max(1,ceil(max(deg)/2));
    x_extract = extractsolution(momentsstructure,options);
end
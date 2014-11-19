function [sol,x_extract,momentsstructure,sosout] = solvemoment(F,obj,options,k)
%SOLVEMOMENT Application of Lasserre's moment-method for polynomial programming
%
%   min        h(x)
%   subject to
%              F(x) >= 0,
%
%   [DIAGNOSTIC,X,MOMENT,SOS] = SOLVEMOMENT(F,h,options,k)
%
%   diagnostic : Struct with diagnostics
%   x          : Extracted global solutions
%   moment     : Structure with various variables needed to recover solution
%   sos        : SOS decomposition {max t s.t h-t = p0+sum(pi*Fi), pi = vi'*Qi*vi}
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

% Get all element-wise constraints, and put them in a vector
% Furthermore, gather the other constraints and place them
% in a new LMI object.
% Additionally, we find out the variables on which we perform
% the relaxation.
vecConstraints = [];
sdpConstraints = [];
isinequality = [];
binaries = [];
xvars = [];
Fnew = ([]);
for i = 1:length(F)
    if is(F(i),'elementwise')
        X = sdpvar(F(i));
        vecConstraints = [vecConstraints;X(:)];
        isinequality = [isinequality ones(1,prod(size(X)))];
        xvars = [xvars depends(X(:))];
    elseif is(F(i),'equality')
        X = sdpvar(F(i));
        if is(X,'symmetric')
            X = X(find(triu(ones(length(X)))));
        end
        vecConstraints = [vecConstraints;-X(:)];
        isinequality = [isinequality zeros(1,prod(size(X)))];
        xvars = [xvars depends(X(:))];
    elseif is(F(i),'sdp')
        sdpConstraints{end+1} = sdpvar(F(i));
        xvars = [xvars depends(F(i))];
    elseif is(F(i),'binary')
        binaries = [binaries getvariables(F(i))];
    else
        Fnew = Fnew+F(i); % Should only be SOCP constraints
    end
end

% Recover the involved variables
x = recover(unique([depends(obj) xvars]));
n = length(x);

% Check degrees of constraints
deg = [];
for i = 1:length(vecConstraints)
    deg(end+1) = degree(vecConstraints(i));
end
for i = 1:length(sdpConstraints)
    deg(end+1) = degree(sdpConstraints{i});
end
if isempty(deg)
    deg = 0;
end

% Perform Schur complements if possible to reduce the
% dimension of the SDP constraints
% for i = 1:length(sdpConstraints)
%     Fi = sdpConstraints{i};
%     j = 1;
%     while j<=length(Fi) & (length(Fi)>1)
%         if isa(Fi(j,j),'double')
%             ks = 1:length(Fi);ks(j)=[];
%             v = Fi(ks,j);
%             vv = v*v'/Fi(j,j);
%             if degree(vv)<=max(deg)
%                 Fi = Fi(ks,ks) - vv;
%             end
%         else
%             j = j+1;
%         end
%     end
% end

% Create lowest possible relaxation if k=[]
% k_min = floor((max(degree(obj)+1,max(deg))+1)/2); % why did I use this?
d = ceil((max(degree(obj),max(deg)))/2);
k_min = d;
if isempty(k)
    k = k_min;
else
    if k<k_min
        error('Higher order relaxation needed')
    end
end

% Generate monomials of order k
u{k} = monolist(x,k);
ulong{k} = monolist(x,2*k);

% Largest moment matrix. NOTE SHIFT M{k+1} = M_k.
M{k+1}=u{k}*u{k}';
% Moment matrices easily generated with this trick
% The matrices will NOT be rank-1 since the products
% generate the relaxed variables

% ... and lower degree localization matrices
M{1} = 1;
for i = 1:1:k-1;
    n_i = round(factorial(n+k-i)/(factorial(n)*factorial(k-i)));
    M{k-i+1} = M{k+1}(1:n_i,1:n_i);
end

% %Moment structure exploition
% M_original = M{end};
% setsdpvar(recover(getvariables(M_original)),0*getvariables(M_original)');
% if 1%options.moment.blockdiag & isempty(F)
%     Ms = blockdiagmoment(obj,M);
% end

% Lasserres relaxation (Lasserre, SIAM J. OPTIM, 11(3) 796-817)
Fmoments = (M{k+1}>=0);
for i = 1:length(vecConstraints)   
    if isinequality(i)
        v_k = floor((degree(vecConstraints(i))+1)/2);
        Localizer = vecConstraints(i)*M{k-v_k+1};
        if isa(vecConstraints(i),'double')
            if vecConstraints(i)<0
                error('Problem is trivially infeasible due to negative constant')
            else
                continue
            end
        end
        Fmoments = Fmoments+(Localizer>=0);
    else
        if isa(vecConstraints(i),'double')
            if vecConstraints(i)~=0
                error('Problem is trivially infeasible due to non-zero constant in equality constraints')
            else
                continue
            end
        end        
        Localizer = vecConstraints(i)*monolist(x,2*k-degree(vecConstraints(i)));      
        Fmoments = Fmoments+(Localizer==0);
    end
end
for i = 1:length(sdpConstraints)
    v_k = floor((degree(sdpConstraints{i})+1)/2);
    Fmoments = Fmoments+(kron(M{k-v_k+1},sdpConstraints{i})>=0);
end

% Add them all
Fnew = Fnew + Fmoments;%unblkdiag(Fmoments);

% No objective, minimize trace on moment-matrix instead
if isempty(obj)
    obj = trace(M{k+1});
end

% Get all binary and reduce problem
binaries = union(binaries,yalmip('binvariables'));
if ~isempty(binaries)
    obj = eliminateBinary(obj,binaries);
    for i = 1:length(Fmoments)
        Fnew(i) = eliminateBinary(Fnew(i),binaries);
    end
    for i = 2:1:k+1;
        M{i} = eliminateBinary(M{i},binaries);
    end
end

% Solve
sol = solvesdp(Fnew,obj,sdpsettings(options,'relax',1));

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
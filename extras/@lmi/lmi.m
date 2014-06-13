function F = lmi(X,handlestring,dummy,noprune,symmetryKnown)
% The command lmi is obsolete, see help sdpvar/set

% Author Johan Löfberg
% $Id: lmi.m,v 1.34 2009-05-29 08:05:12 joloef Exp $

superiorto('sdpvar')

fast = 0;
noprune = 0;
symmetryKnown = 0;
switch nargin
    case 0
        F.clauses={};
        F.LMIid    = [];
        F.savedata = [];
        F = class(F,'lmi');
        return
    case 1
        handlestring = '';
        % ok probably, check later
    case 2
        if ~isa(handlestring,'char')
            error('The handle must be a string')
        end
    case 4
        handlestring = '';
        noprune = 1;
    case 5
        handlestring = '';
        symmetryKnown = 1;

    otherwise
        error('Wrong number of arguments')
end

F.clauses  = {};
F.LMIid    = [];
F.savedata = [];

try
    F = class(F,'lmi');
catch
    disp('Failure in creating SET object. Restart MATLAB and try again, or type clear classes')
end

% Return empty LMI in degenerate case
if isempty(X)
    % A bit inconsitently, we create a constraint, which is empty. We count
    % the number of ocnstraints by looking at LMIid, which still is zero
    F.clauses{1}.data=[];
    F.clauses{1}.type = [];
    F.clauses{1}.symbolic=[];
    F.clauses{1}.handle=[];
    F.clauses{1}.strict = [];
    F.clauses{1}.cut = [];
    F.clauses{1}.expanded = [];
    F.clauses{1}.lift = [];    
    F.clauses{1}.schurfun  = [];
    F.clauses{1}.schurdata = [];
    F.clauses{1}.jointprobabilistic = [];
    F.clauses{1}.confidencelevel = [];
    F.clauses{1}.extra = [];
    return
end

% 0 : Nothing (i.e error)
% 1 : LMI
% 2 : Elementwise
% 3 : Equality
% 4 : SOCC
% 5 : Rotated Lorentz
% 9 : KYP constraint
%10 : Eigenvalue constraint
%11 : Sum of square constraint
%12 : Logic CNF
%15 : Deterministic uncertain variable
%16 : Random uncertain variable
%20 : Power cone
%30 : User generated Schur
%40 : Generalized KYP
%50 : SOS2
%51 : SOS1
%52 : semivar
%53 : semiintvar
%55 : Complementary
%56 : Meta constraint to be expanded (implies, iff)
%60 : Chance constraint

switch class(X)
    case 'lmi'
        F = X;
        return
    case 'char'
        try
            % Parse (check for old notation, take care of ||...|| notation)
            sdpvarExpr = parseLMI(X);
            % Old notation
            X = strrep(X,'.>','>');
            X = strrep(X,'.<','<');            
        catch
            error(lasterr)
        end
        % Evaluate
        try
            Fi=evalin('caller',sdpvarExpr);
            switch class(Fi)
                case 'sdpvar'
                    Fi = {Fi};
                    strict = 0;
                case 'constraint'
                    [Fi,strict] = getlist(Fi);
                case 'lmi'
                    F = Fi;
                    return;                    
                otherwise
                    error('The string does not define an SDPVAR object')
            end

        catch
            error(lasterr)
        end

    case 'sdpvar'
        Fi = {X};
        strict = 0;
        X='Numeric value';

    case 'logic'
        Fi = {X};
        strict = 0;
        X = 'Numeric value';

    case {'struct','constraint'}      
        [Fi,strict,LMIIdentifiers,tags] = getlist(X); 
        if isempty(handlestring) || length(handlestring)==0
            handlestring = tags{1};
        end
        X='Numeric value';

    otherwise
        error('The first argument should be an sdpvar object or string.')
end

% Return empty LMI in degenerate case
if all(cellfun('isempty',Fi))
    return
end

TypeofConstraint = zeros(length(Fi),1)-1;

% support for SOS constraints placed in vector
if length(Fi) == 1
    if isa(Fi{1},'sdpvar')% is(Fi{1},'sos')
        if gethackflag(Fi{1})==1
            % Expand to a set of SOS constraints
            if ~issymmetric(Fi{1})
                p = Fi{1}(:);
                for i = 1:length(p)
                    Fi{i} = p(i);
                end
                TypeofConstraint = zeros(length(Fi),1)-1;
                strict = zeros(length(Fi),1)-1;
            end
        end
    end
end

i = 1;
while i <= length(Fi)
    thisFi = Fi{i};
    if ~isempty(thisFi)
        if isa(thisFi,'logic')
            TypeofConstraint(i) = 12;
        else
            TypeofConstraint(i) = gethackflag(Fi{i});
        end
        % User just sent a sdpvar object, interpret as X>0
        if (TypeofConstraint(i)==0)
            TypeofConstraint(i) = 1;
        end

        if (TypeofConstraint(i)==1) && (symmetryKnown == 0)        
            [n,m]=size(thisFi);
            if (n~=m) || ((TypeofConstraint(i) == 1) && (n*m==1)) || ~ishermitian(thisFi)
                TypeofConstraint(i) = 2;
            end
        end

        if TypeofConstraint(i) == 2

            % remove constraint of the type set(0 >= 0)
            B = getbase(thisFi);            
            if ~noprune
                Bv = B(:,2:end);
                notused = find((~any(Bv,2)) & (B(:,1)>=0));
                if ~isempty(notused)
                    used = setdiff(1:size(Bv,1),notused);                    
                    thisFi = thisFi(used);                    
                end
            end
        end

        % save cleaned version
        Fi{i} = thisFi;

        switch TypeofConstraint(i)
            case {1,2,3,4,5,7,8,9,10,11,12,13,15,16,20,30,40,50,51,52,53,54}
                i = i + 1;
            otherwise
                error('Error in argument in LMI. Please report bug');
        end
    end
end

% Special case for x<y<z<...
% Merge into one constraint, gives nicer displays
if ~exist('LMIIdentifiers','var')
   LMIIdentifiers = yalmip('lmiid');
end

if all(TypeofConstraint == 2) && all(strict==strict(1))
    if length(Fi)>1
        vecF = [];
        sizes = zeros(length(Fi),1);
        % Speed up concatenation of merging
        for i = 1:length(Fi)
            sizes(i) = prod(size(Fi{i}));
        end
        if all(sizes == 1)
            vecF = [Fi{:}]';
        else
            for i = 1:length(Fi)
                fi = Fi{i};
                if sizes(i) > 1                    
                    fi = reshape(fi,prod(size(fi)),1);
                end
                vecF = [vecF;fi];
            end
        end
    else
        vecF = reshape(Fi{1},prod(size(Fi{1})),1);
    end
    F.clauses{1}.data=vecF;
    F.clauses{1}.type = 2;
    F.clauses{1}.symbolic=X;
    F.clauses{1}.handle=handlestring;
    F.clauses{1}.strict = strict(1);
    F.clauses{1}.cut = 0;
    F.clauses{1}.expanded = 0;
    F.clauses{1}.lift = 0;
    F.clauses{1}.schurfun  = '';
    F.clauses{1}.schurdata = [];
    F.clauses{1}.jointprobabilistic = [];
    F.clauses{1}.confidencelevel = [];
    F.clauses{1}.extra = [];
    F.LMIid = [F.LMIid LMIIdentifiers(1)];
else
    for i = 1:length(Fi)
        switch TypeofConstraint(i)
            case {1,2,3,4,5,7,8,9,10,11,12,13,15,16,20,30,40,50,51,52,53,54}
                F.clauses{i}.data=Fi{i};
                F.clauses{i}.type = TypeofConstraint(i);
                F.clauses{i}.symbolic=X;
                F.clauses{i}.handle=handlestring;
                F.clauses{i}.strict = strict(i);
                F.clauses{i}.cut = 0;
                F.clauses{i}.expanded = 0;
                F.clauses{i}.lift = 0; 
                F.clauses{i}.schurfun  = '';
                F.clauses{i}.schurdata = [];
                F.clauses{i}.jointprobabilistic = [];
                F.clauses{i}.confidencelevel = [];
                F.clauses{i}.extra = [];
                F.LMIid = [F.LMIid LMIIdentifiers(i)];
                i = i + 1;
            otherwise
                error('Error in argument in LMI. Please report bug');
        end
    end
end


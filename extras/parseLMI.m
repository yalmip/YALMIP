function sdpvarExpr = parseLMI(X)
%PARSELMI Internal function for dirty parsing of lmi string

% Author Johan Löfberg
% $Id: parseLMI.m,v 1.3 2004-07-27 14:11:57 johanl Exp $


% Assume Error

% TypeofConstraint = 1;

% Check for obsolete notation .<, .>, =
single_equality = findstr(X,'=');
if isempty(findstr(X,'>=')) & isempty(findstr(X,'<=')) & (rem(length(single_equality),2) == 1) % There is a single =    
    error('Obsolete constraint =, use == in equalities.');
    %X = strrep(X,'=','==');
    %X = strrep(X,'====','=='); % Whoops, == replaced with =====!
end

if ~isempty(findstr(X,'.<'))
    disp(' ');
    disp('Warning: Obsolete constraint .<, use < in equalities.');
    disp('If you have a Hermitian matrix which you want to')
    disp('constrain to have negative values, use ,e.g.,  P(:)<0')
    disp('The constraint has been changed to <')
    disp(' ');
    X = strrep(X,'.<','<');
end

if ~isempty(findstr(X,'.>'))
    disp(' ');
    disp('Warning: Obsolete constraint .>, use > in equalities.');
    disp('If you have a Hermitian matrix which you want to')
    disp('constrain to have negative values, use e.g.  P(:)>0')
    disp('The constraint has been changed to >')
    disp(' ');
    X = strrep(X,'.>','>');
end

% Any norm? If not, we're done!
if isempty(findstr(X,'||'))
    sdpvarExpr = X;
    return
end


% The only actual parsing needed is for the ||Axplusb||<cx+d notation
delimiters = {'.<','.>','<.','>.','<','>','==',''};
delindex = 0;indtodel = [];
while isempty(indtodel) & (delindex<=7)
    delindex = delindex+1;
    indtodel = findstr(X,delimiters{delindex});
end

switch delindex
    case 5 %<
        TypeofConstraint = 1;
        ind_start =  indtodel-1;
        ind_end   =  indtodel+1;
        REVERSE   = 1;
    case 6 %>
        TypeofConstraint = 1;
        ind_start =  indtodel-1;
        ind_end   =  indtodel+1;
        REVERSE   = 0;
    case 7%=
        TypeofConstraint = 3;
        ind_start = indtodel-1;
        ind_end   = indtodel+2;
        REVERSE   = 0;
    case {3,4}
        error('Incorrect argument. Perhaps you mean .> or .<');
    otherwise
        error('Incorrect argument. Could not find <=>.>.<')
end

LeftHand   = X(1:ind_start);
RightHand  = X(ind_end:end);

if REVERSE
    temp = LeftHand;
    LeftHand = RightHand;
    RightHand = temp;
end

% Search for a norm expression
ind_norm_Right = findstr(RightHand,'||');
ind_norm_Left  = findstr(LeftHand,'||');

% Any norm at all, if not, we're done!
if isempty(ind_norm_Right) & isempty(ind_norm_Left)
    sdpvarExpr = [LeftHand '-(' RightHand ')'];
    return
end

% Equality constrained norm?
if (TypeofConstraint == 3) & (ind_norm_Right | ind_norm_Left)
    error('Equality constraints cannot be used with ||..|| operator');
end

% Convex normconstraint?
if ~isempty(ind_norm_Left)
    error('Norm ||..|| must look like || column vector ||< scalar, or scalar>|| vector ||')
end

% Everthing seem ok in norm
TypeofConstraint = 4;
Axplusb = RightHand(2+ind_norm_Right(1):ind_norm_Right(2)-1);
WithoutNorm = strrep(RightHand,RightHand(ind_norm_Right(1):ind_norm_Right(2)+1),'0');
if length(WithoutNorm)==1
    cxplusd = LeftHand;
else
    cxplusd = [LeftHand '-(' WithoutNorm ')'];
end

sdpvarExpr = ['cone( ' Axplusb ',' cxplusd ')'];



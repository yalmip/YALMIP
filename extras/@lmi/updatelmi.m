function sys = updatelmi(F,thelmi)
%updatelmi         (OBSOLETE) Updates the numerical content of a symbolic LMI
%   
%    F = UPDATELMI(F,n)         Updates the n'th symbolic LMIs
%    F = UPDATELMI(F,'tag')     Updates the symbolic LMIs having specified tag
%
%    See also   LMI, DELLMI


% Author Johan Löfberg
% $Id: updatelmi.m,v 1.2 2004-07-19 13:54:37 johanl Exp $

error('updatelmi is obsolete');

if nargin ~=2
    error('UPPDATELMI needs two argument')
end

if ~(isa(F,'lmi') & (isa(thelmi,'double') | isa(thelmi,'char')))
    error('First argument should be an lmi object and second argument integer or string')
end

% If indexed using handle, convert to index
if isa(thelmi,'char')
    thelmitemp = [];
    for i = 1:size(F.AllConstraints,2)
        if strcmp(F.AllConstraints{i}.handle,thelmi)
            thelmitemp=[thelmitemp i];
        end
    end
    %     for i = 1:size(F.Equalities,2)
    %       if strcmp(F.Equalities{i}.handle,thelmi)
    % 	thelmitemp=[thelmitemp i+size(F.AllConstraints,2)];
    %       end
    %     end	
    if isempty(thelmitemp)
        em = ['LMI ''' thelmi ''' not available.'];
        error(em)
    else
        thelmi = thelmitemp;
    end
end

% Checks so that it exist
if any((thelmi<1)) | any((thelmi>size(F.AllConstraints,2)))
    em = ['LMI #' num2str(thelmi) ' not available.'];
    error(em)
end

for j = 1:length(thelmi)
    nlmi = thelmi(j);
    % Get the old symbolic expression
    if (nlmi>size(F.AllConstraints,2))
        % An equality constraint 
        nlmi = nlmi-size(F.AllConstraints,2);
        X = F.Equalities{nlmi}.symbolic;
        eq = 1;
    else
        % LMI
        X = F.AllConstraints{nlmi}.symbolic; 
        eq = 0;
    end
    
    % Check for symbolic value
    if strcmpi('Numeric value',X)
        error('Numeric LMIs cannot be updated');
    end
    
    % Parse the LMI
    sdpvarExpr =  parseLMI(X); 
    X = strrep(X,'.>','>');
    X = strrep(X,'.<','<'); 
    X = strrep(X,'=','==');X = strrep(X,'====','==');    
    try
        Fi=evalin('caller',sdpvarExpr);
    catch
        error(lasterr)
    end
    
    TypeofConstraint = gethackflag(Fi);  
    % User just sent a sdpvar object, interpret as X>0
    if (TypeofConstraint==0)
        TypeofConstraint = 1;
    end
    
    % check so that LMIs are symmetric
    if (TypeofConstraint == 1)
        [n,m]=size(Fi);
        if (n~=m) | (~issymmetric(Fi)) |(n*m==1)
            TypeofConstraint = 2; %Change to elementwise	
        end
    end
    
    % Check the normbound
    if (TypeofConstraint == 4)
        [n,m]=size(Fi);
        if m~=1
            error('Norm ||..|| must look like || vector ||< scalar, or scalar>|| vector ||')	
        end
    end
    
    sys=F;
    switch TypeofConstraint
        case {1,2,3,4}
            sys.AllConstraints{nlmi}.data=Fi;
            sys.AllConstraints{nlmi}.type = TypeofConstraint;
            sys.AllConstraints{nlmi}.symbolic=X;
        otherwise
            error('Error in argument.');
    end
end;




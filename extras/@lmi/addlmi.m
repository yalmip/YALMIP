function sys = addlmi(F,X,handlestring)
% ADDLMI Obsolete. 
%
% Use overloaded + instead.
%

% Author Johan Löfberg
% $Id: addlmi.m,v 1.2 2004-07-19 13:54:35 johanl Exp $


error('addlmi is obsolete. Use overloaded +');

switch(nargin)
    case {0,1}
        error('ADDLMI needs at least two arguments')
    case 2
        %       sys = F + lmi(X);
        handlestring = '';
    case 3
        %       sys = F + lmi(X,handlestring);
        if ~isa(handlestring,'char')
            error('The handle must be a string')
        end
    otherwise
        error('addlmi takes at mopst 3 arguments');
end

% Return empty LMI in degenerate case
if isempty(X)
    sys = F;
    return
end

%Symbolic expression given 
% 0 : Nothing (i.e. error)
% 1 : LMI
% 2 : Elementwise
% 3 : Equality
TypeofConstraint = 0;
switch class(X)
    case 'char'
        try
            sdpvarExpr =  parseLMI(X);
            X = strrep(X,'.>','>');
            X = strrep(X,'.<','<'); 
            X = strrep(X,'=','==');X = strrep(X,'====','==');
        catch
            error(lasterr)
        end
        try
            Fi=evalin('caller',sdpvarExpr); 
            if ~(isa(Fi,'sdpvar') | isa(Fi,'constraint'))
                error('The string does not define a constraint involving SDPVAR objects')
            end
        catch
            error(lasterr)
        end
    case 'sdpvar'
        Fi = X;
        X='Numeric value';
    case 'constraint'
        Fi = getlist(X);
        if length(Fi)>1
            error('Hey, you have updated your code for version 3, so why don''t you remove those addlmi''s');
        else
            Fi = Fi{1};
        end
        X='Numeric value';       
    otherwise
        error('LMI should be an sdpvar object or string.')
end

if isempty(Fi)
    sys = F;
    return
end

switch nargin
    case 2
        sys = F + lmi(Fi);
    case 3
        sys = F + lmi(Fi,handlestring);
    otherwise
end
return
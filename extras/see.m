function see(X,showfull)
%SEE Displays internal structure of variable
%
%    Shows the base matrices that build up the variable
%  
%    X = X0+x_1*X1+x_2*X2+...
%  
%    SEE(X)            Display matrices Xi in dense format
%    SEE(X,'sparse')   Display matrices Xi in sparse format
%
%    See also   SDPVAR

disp(' ');
if nargin==1  
    disp('Constant matrix');disp(' ')
    disp(full(X))
    disp('Base matrices');disp(' ')
    disp('[]')
    disp(' ')			
    disp('Used variables');disp(' ')
    disp('[]')
else	
    switch showfull
        case 'sparse'
            disp('Constant matrix');disp(' ')
            disp((X))
            disp('Base matrices');disp(' ')
            disp('[]')
            disp(' ')			
            disp('Used variables');disp(' ')
            disp('[]')
            
        case 'full'
            see(X)
        otherwise
            error('Second argument should be ''sparse'' or ''full')
    end
end
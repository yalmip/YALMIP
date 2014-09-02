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
    disp(full(getbasematrix(X,0)))
    disp('Base matrices');disp(' ')
    for i = 1:length(X.lmi_variables);       
        disp(full(getbasematrix(X,X.lmi_variables(i))))       
        disp(' ')			
    end;
    disp('Used variables');disp(' ')
    disp(X.lmi_variables)
else	
    switch showfull
        
        case 'sparse'
            disp('Constant matrix');disp(' ')
            disp((getbasematrix(X,0)))
            disp('Base matrices');disp(' ')
            for i = 1:length(X.lmi_variables);
                disp((getbasematrix(X,X.lmi_variables(i))));
                disp(' ')
            end;
            disp('Used variables');disp(' ')
            disp(X.lmi_variables)
            
        case 'full'
            see(X)
            
        otherwise
            error('Second argument should be ''sparse'' or ''full')
            
    end
end
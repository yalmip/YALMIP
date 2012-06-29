function define(varargin)
%DEFINE
% 
% definestring = 'sdpvar ';
% 
% for i=1:nargin
%     X = varargin{i};
%     xname = inputname(i);
%     [n,m]=size(X);        
%     for ii = 1:n
%         for jj=1:m
%             x=X(ii,jj);
%             v=getvariables(x);
%             if min(n,m)>1
%                 dostring = [xname num2str(ii) num2str(jj) '=recover(' num2str(v) ');'];
%             else
%                 dostring = [xname num2str(max(ii,jj)) '=recover(' num2str(v) ');'];
%             end
%             evalin('caller',dostring);
%         end
%     end       
% end

for i=1:nargin
    X = varargin{i};
    xname = inputname(i);
    [n,m]=size(X);   
    namesout = [];
    variablesin = [];
    for ii = 1:n
        for jj=1:m
            x=X(ii,jj);
            v=getvariables(x);
            if min(n,m)>1
                namesout = [namesout xname num2str(ii) num2str(jj) ','];
                variablesin = [variablesin num2str(v) ' '];
             %   dostring = [xname num2str(ii) num2str(jj) '=recover(' num2str(v) ');'];
            else
                variablesin = [variablesin num2str(v) ' '];
                namesout = [namesout xname num2str(max(ii,jj)) ','];
             %   dostring = [xname num2str(max(ii,jj)) '=recover(' num2str(v) ');'];
            end
           % evalin('caller',dostring);
        end
    end          
    dostring = ['[' namesout(1:end-1) '] = recover([ ' variablesin ']);'];
    evalin('caller',dostring);    
end
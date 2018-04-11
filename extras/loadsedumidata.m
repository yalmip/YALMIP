function [F,h] = loadsedumidata(varargin)
%LOADSEDUMIDATA Loads a problem definition in the SeDuMi format
%
%    [F,h] = loadsedumidata('filename')  Loads the problem min(h(x)), F(x)>0 from file 'filename'
%    [F,h] = loadsedumidata              An "Open" - box will be opened

filename = varargin{1};

% Does the file exist
if ~exist(filename)
    filename = [filename '.mat'];
    if ~exist(filename)
        error(['No such file.']);
    end
end
   
load(filename)
try   
    if ~exist('At')
        At = A;
    end
    if ~exist('b')
        b = zeros(size(At,1),1);
    else
        b = b(:);
    end
    if ~exist('c')
        if exist('C')
            c = C(:);
        else
            c = zeros(size(At,2),1);
        end
    else
        c = c(:);
    end
    
    K = K;
catch
    error('The file should contain the data At, b, c and K');
end

[F,h] = sedumi2yalmip(At,b,c,K);
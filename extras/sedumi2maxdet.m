function [F_struc,F_blksz,G_struc,G_blksz]  = sedumi2maxdet(F_struc,K)
%SEDUMI2MAXDET Internal function to convert SeDuMi structure to format needed in MAXDET

% Author Johan Löfberg
% $Id: sedumi2maxdet.m,v 1.4 2006-12-18 14:42:28 joloef Exp $

switch K.m(1)
    case 0
        % No MAXDET terms
        G_struc = [];
        G_blksz = [];
        F_struc = F_struc;
        if any(K.s>0)
            F_blksz = [repmat(1,1,K.l) K.s];
        else
            F_blksz = [repmat(1,1,K.l)];
        end
    case 1
        % Error, FIXME
        G_struc = F_struc(K.l,:);
        G_blksz = [1];
        F_blksz = [repmat(1,1,K.l-1) K.s];
        F_struc = [F_struc(1:1:K.l-1,:);F_struc(K.l+1:1:end,:)];
    otherwise
        % Number of maxdet terms
        m = length(K.m);
        % #rows for these matrixes
        r = sum((K.m).^2);
        G_struc = F_struc(end-r+1:end,:);
        G_blksz = K.m;
        if length(K.s) > length(K.m)
            F_blksz = [repmat(1,1,K.l) K.s(1:end-m)];
        else
            F_blksz = repmat(1,1,K.l);
        end
        F_struc = F_struc(1:end-r,:);
end

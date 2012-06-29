function monom = recovermonoms(newton_m,x)
%RECOVERMONOMS Internal function used in SOS programs

% Author Johan Löfberg
% $Id: recovermonoms.m,v 1.1 2006-08-10 18:00:22 joloef Exp $

if size(newton_m,1)==1 & nnz(newton_m)==0
    monom = 1;
    return
end

[mt,oldvariabletype,mt_hash,hash] = yalmip('monomtable');

origSize = size(mt,1);
vars = x.lmi_variables;

newton_m_extended = spalloc(size(newton_m,1),size(mt,2),nnz(newton_m));
newton_m_extended(:,vars) = newton_m;
newton_m_extended_hash = newton_m_extended*hash;

newton_m = newton_m';
newton_m_extended = newton_m_extended';
mt = mt';

monom = x;
%monom.n = size(newton_m,1);
monom.dim(1) = size(newton_m,2);
monom.dim(2) = 1;
monom.lmi_variables = [];
%variable_here = ones(1,size(newton_m,1));
%variable_here = ones(1,size(newton_m,2));
mt_hash = full(mt_hash);
nvar = size(mt,2);
new_mt_hash = [];
new_mt = [];
variable_here = any(newton_m,1);
for i = find(variable_here)%1:size(newton_m,2)
    %if nnz(newton_m(:,i))==0
    %    variable_here(i) = 0;
    %else
        if isempty(mt_hash)%MESSY DUE TO BEHAVIOUR IN LINUX 6.1
           previous_variable = [];
        else
            previous_variable = find(mt_hash == newton_m_extended_hash(i));            
        end
        if isempty(previous_variable)  
            if isempty(new_mt_hash)
                previous_variable = [];
            else
                previous_variable = find(new_mt_hash == newton_m_extended_hash(i));
            end
            if isempty(previous_variable)
                nvar = nvar + 1;
          %      mt = [mt newton_m_extended(:,i)];
          %      mt_hash = [mt_hash ; newton_m_extended(:,i)'*hash];
                new_mt = [new_mt newton_m_extended(:,i)];
                new_mt_hash = [new_mt_hash ; newton_m_extended(:,i)'*hash];
          
                monom.lmi_variables = [monom.lmi_variables nvar];
            else
                monom.lmi_variables = [monom.lmi_variables size(mt,2)+previous_variable];
            end
        else
%            try
                 monom.lmi_variables = [monom.lmi_variables previous_variable];
 %           catch
  %              error
   %         end
        end
  %  end
end
mt = [mt new_mt]';
mt_hash = [mt_hash;new_mt_hash];
%mt = mt';

% if append_one    
%     monom.basis = [[1;spalloc(monom.dim(1),1,0)] [spalloc(1,monom.dim(1),0);speye(monom.dim(1))]];
%     monom.dim(1) = monom.dim(1) + 1;
% else
% monom.basis = [spalloc(monom.dim(1),1,0) spalloc(monom.dim(1),length(find(variable_here)),0)];
% monom.basis = [spalloc(monom.dim(1),1+length(find(variable_here)),1+length(find(variable_here)))];
% monom.basis(find(variable_here),2:end) = speye(length(find(variable_here)));
% monom.basis(find(variable_here==0),1)=1;

nz = find(variable_here);
zv = find(variable_here==0);
i = [zv(:);nz(:)];
j = [ones(length(zv),1);(2:(length(nz)+1))'];
k = ones(length(variable_here),1);
monom.basis = sparse(i,j,k,monom.dim(1),1+length(find(variable_here)));

%end

% Fucked up order
if any(diff(monom.lmi_variables)<0)
    [i,j]=sort(monom.lmi_variables);
    monom.lmi_variables = monom.lmi_variables(j);
    monom.basis(:,2:end) = monom.basis(:,j+1);
end

un_monom_vars = uniquestripped(monom.lmi_variables);
if length(un_monom_vars)<length(monom.lmi_variables)
    [un_monom_vars,hh,jj] = unique(monom.lmi_variables);
    if length(monom.lmi_variables) ~=length(un_monom_vars)
        newmonombase = monom.basis*sparse([1 1+jj],[1 1+(1:length(jj))],ones(1,1+length(jj)))';
        monom.basis = newmonombase;
        monom.lmi_variables = un_monom_vars;
    end
end

if size(mt,1) > origSize
    newmt = mt(origSize+1:end,:);
    newvariabletype = spalloc(size(newmt,1),1,0)';
    nonlinear = ~(sum(newmt,2)==1 & sum(newmt~=0,2)==1);
    if ~isempty(nonlinear)
        %mt = internal_sdpvarstate.monomtable;
        newvariabletype(nonlinear) = 3;
        quadratic = sum(newmt,2)==2;
        newvariabletype(quadratic) = 2;
        bilinear = max(newmt,[],2)<=1;
        newvariabletype(bilinear & quadratic) = 1;
        sigmonial = any(0>newmt,2) | any(newmt-fix(newmt),2);
        newvariabletype(sigmonial) = 4;
    end
    yalmip('setmonomtable',mt,[oldvariabletype newvariabletype]);
end



function sys = pretty(F)
% PRETTY Pretty print an SET object
%    PRETTY(F) tries to print the SET F in a format that
%    resembles type-set mathematics.
% 
% Requires the Symbolic Math Toolbox
  
% Author Johan Löfberg
% $Id: pretty.m,v 1.3 2005-02-04 10:10:27 johanl Exp $
   
error('Pretty is obsolte and no longer supported')

  nlmi = size(F.clauses,2);
  
  spaces = ['                                    '];
  if (nlmi == 0) & (neq == 0)
    disp('empty LMI')
    return
  end
  
  lmiinfo{1} = 'Matrix inequality';
  lmiinfo{2} = 'Element-wise inequality';
  lmiinfo{3} = 'Equality constraint';
  lmiinfo{4} = 'Second order cone constraint';
  
  % Needed for formating
  infolen    = max(cellfun('length',lmiinfo));
  
  
  % Check if there are any handles, sizes of strings
  max_ns = 0;
  max_nh = 0;
  for i = 1:nlmi
    max_ns = max(max_ns,length(F.clauses{i}.symbolic));
    max_nh = max(max_nh,length(F.clauses{i}.handle));
  end
  
  ns = min(40,max_ns);
  nh = min(40,max_nh);
  nt = infolen+3;
  nh = 40;
  
  if nlmi>0
    nn = ceil(log10(nlmi+1)+1); 
    disp('*************************');
    disp('Inequality constraints');
    disp('*************************'); 
    for i = 1:nlmi
	type     = lmiinfo{F.clauses{i}.type};
	handle   = truncstring(F.clauses{i}.handle,nh);
	number   = strjust(fillstring(['#' num2str(i)],nn),'right');
	if isempty(handle)
	  infostr = sprintf('%s   %s',number,type);
	else
	  infostr = sprintf('%s   %s    (%s)',number,type,handle);
	end
	disp(infostr);
	prettyhack(F.clauses{i}.symbolic);
	disp(' ' );
    end
    disp(' ');
  else
    disp('*************************');  
    disp('No inequality constraints');
    disp('*************************');
  end
  
  
function x= truncstring(x,n)
  if length(x) > n
    x = [x(1:n-3) '...'];
  end
  
function x = fillstring(x,n)
  x = [x blanks(n-length(x))];
  
  % Works for some type of problems
function prettyhack(F)
  
% Some hacks to make it work
  Y = F;
  F = strrep(F,'''','^T');
  F = strrep(F,'diag','Diag');
  F = strrep(F,'min','Min');
  F = strrep(F,'max','Max');
  
  try
    TypeofConstraint = 0;  
    delimiters = {'.<','.>','<.','>.','<','>','=',''};  
    delindex = 0;indtodel = [];
    while isempty(indtodel) & (delindex<=7)
      delindex = delindex+1;
      indtodel = findstr(F,delimiters{delindex});   
    end
    switch delindex
     case 5 %<
      TypeofConstraint = 1;
      ind_start =  indtodel-1;
      ind_end   =  indtodel+1;
      REVERSE   = 1;
      dell = '<';
     case 6 %>
      TypeofConstraint = 1;
      ind_start =  indtodel-1;
      ind_end   =  indtodel+1;
      REVERSE   = 0;
      dell = '>';
     case 1 %.<
      TypeofConstraint = 2;
      ind_start =  indtodel-1;
      ind_end   =  indtodel+2;
      REVERSE   = 1;
      dell = '<';
     case 2 %.>
      TypeofConstraint = 2;
      ind_start =  indtodel-1;
      ind_end   =  indtodel+2;
      REVERSE   = 0;
      dell = '>';
     case 7%=
      TypeofConstraint = 3;
      ind_start = indtodel-1;
      ind_end   = indtodel+1;
      REVERSE   = 0;
      dell = '=';
     case {3,4}
      error('Incorrect argument. Perhaps you mean .> or .<');   
     otherwise
      error('Incorrect argument. Could not find <=>.>.<')
    end
    
    LeftHand   = char(sym(F(1:ind_start)));
    RightHand  = char(sym(F(ind_end:end)));
    
    maplemex('interface(prettyprint=true);');
    disp(' ')
    r = maplemex(['print(' LeftHand ' ' dell ' ' RightHand ');']);
    if ~isempty(r)
      disp(' ')
      disp([strjust(fillstring(Y,79),'center') ' (cannot be displayed in Maple)'])
    end
  catch
    disp(' ');
    disp([strjust(fillstring(Y,79),'center') ' (cannot be displayed in Maple)'])
  end 
  maplemex('interface(prettyprint=false);');
  
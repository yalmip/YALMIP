function YESNO = is(X,property,additional)
%IS Check property of variable.
%   d = IS(x,property) returns 1 if 'property' holds
%
%   Properties possible to test are: 'real', 'symmetric', 'hermitian',
%   'scalar', 'linear', 'bilinear','quadratic','sigmonial', 'homogeneous', 'integer', 'binary'

switch property
    case 'logic'
        YESNO = X.typeflag==12;
    case 'binary'
        YESNO = any(ismember(depends(X),yalmip('binvariables')));
    case 'integer'
        YESNO = any(ismember(depends(X),yalmip('intvariables')));
    case 'quantized'              
        YESNO = all(ismember(depends(X),yalmip('quantvariables')));
    case 'real'
        YESNO = isreal(X);
    case 'complex'
        YESNO = ~isreal(X.basis);
    case 'interval'
        YESNO = isa(X.basis,'intval');        
    case 'symmetric'
        YESNO = issymmetric(X);
    case 'hermitian'
        YESNO = ishermitian(X);
    case 'scalar'
        YESNO = prod(X.dim)==1;
    case 'compound'
        YESNO = any(ismember(getvariables(X),yalmip('extvariables')));       
    case 'linear'
        variabletype = yalmip('variabletype');
        variabletype = variabletype(X.lmi_variables);       
        YESNO = ~any(variabletype);              
    case 'bilinear'
        variabletype = yalmip('variabletype');
        variabletype = variabletype(X.lmi_variables);                
        YESNO  = all(variabletype<=1);
     case 'quadratic'
         variabletype = yalmip('variabletype');
         variabletype = variabletype(X.lmi_variables);
         YESNO = all(variabletype<=2);
         
    case 'LBQS'
        % Fast code for use in display etc.
        % Checks linearity, bilinearity etc in one call.
        quadratic = 0;
        bilinear  = 0;
        linear    = 0;
        sigmonial = 0;
              
        variabletype = yalmip('variabletype');
        variabletype = variabletype(X.lmi_variables);
        
        linear    = all(variabletype==0);        
        bilinear  = all(variabletype<=1) && ~linear;        
        quadratic = all(variabletype<=2) && ~bilinear;
        sigmonial = any(variabletype==4);
        
        YESNO = full([linear bilinear quadratic sigmonial]);
         
         
         
    case 'lpcone'
          base = X.basis;
          YESNO = all(base(:,1)==0); % No constant
          base = base(:,2:end);
          YESNO = YESNO && all(sum(base,2)==1); 
          YESNO = YESNO && all(sum(base,1)==1);
          YESNO = YESNO && all(sum(base~=0,2)==1);
          YESNO = YESNO && all(sum(base~=0,1)==1);
          YESNO = full(YESNO);

    case 'shiftlpcone'
          base = X.basis;          
          base = base(:,2:end);
          YESNO = all(sum(base,2)==1); 
          YESNO = YESNO && all(sum(base,1)==1);
          YESNO = YESNO && all(sum(base~=0,2)==1);
          YESNO = YESNO && all(sum(base~=0,1)==1);
          YESNO = full(YESNO);

    case 'complexsdpcone'
        if isequal(X.conicinfo,[sqrt(-1) 0])
            YESNO = 1;
            return
        else
            YESNO = 0;
        end
        
    case 'realsdpcone'
        if isequal(X.conicinfo,[1 0])
            YESNO = 1;
            return
        else
            YESNO = 0;
        end
        
    case 'sdpcone'
        if isequal(X.conicinfo,[1 0]) || isequal(X.conicinfo,[sqrt(-1) 0])
            YESNO = 1;
            return
        end
        base = X.basis;
        n = X.dim(1);
        YESNO = full(issymmetric(X) && nnz(base)==n*n && all(sum(base,2)==1) && all(base(:,1)==0)) && length(X.lmi_variables)==n*(n+1)/2 && isreal(base); 
        
        % Fixes bug #15
        if length(X.lmi_variables)>1
            if ~all(diff(X.lmi_variables)==1)
                YESNO=0;
                return;
            end
        end
                
    case 'shiftsdpcone'
        
        if isequal(X.conicinfo,[1 0])
            YESNO = 1;
            return
        elseif isequal(X.conicinfo,[1 1])
            YESNO = 1;
            return
        elseif  isequal(X.conicinfo,[sqrt(-1) 0])
            YESNO = 1;
            return
        elseif isequal(X.conicinfo,[1 1])
            YESNO = 1;
            return
        end
        
        % Fixes bug #15
        if length(X.lmi_variables)>1
            if ~all(diff(X.lmi_variables)==1)
                YESNO=0;
                return;
            end
        end
        
        base = X.basis;
        n = X.dim(1);
        base(:,1)=0;
        YESNO = full(issymmetric(X) && nnz(base)==n*n && all(sum(base,2)==1)) && length(X.lmi_variables)==n*(n+1)/2 && isreal(X);
        if YESNO
            % Possible case
            % FIX : Stupidly slow and complex
            [i,j,k] = find(base');
            Y = reshape(1:n^2,n,n);
            Y = tril(Y);
            Y = (Y+Y')-diag(sparse(diag(Y)));
            [uu,oo,pp] = unique(Y(:));
            YESNO = isequal(i,pp+1);            
        end
        
    case 'socone'
        base = X.basis;
        n = X.dim(1);
        YESNO = X.dim(1)>1 && X.dim(2)==1 && length(X.lmi_variables)==n;
        if YESNO
            cb = base(:,1);
            vb = base(:,2:end);            
            YESNO = YESNO && (nnz(cb)==0) && (nnz(vb-speye(n))==0);
        end   
        
    case 'sigmonial'
          monomtable = yalmip('monomtable');
          monomtable = monomtable(getvariables(X),:);
          YESNO = any(find(any(0>monomtable,2) | any(monomtable-fix(monomtable),2)));   
          
    case 'general'
        evalvariables = yalmip('extvariables');
        YESNO = ~isempty(intersect(getvariables(X),evalvariables));
        
    case 'nonlinear'
        YESNO = ~islinear(X);   
        
    case 'homogeneous'    
        [sqrList,CompressedList] = yalmip('nonlinearvariables');
        [LinearTerms,NonlinearVariables] = getvariables(X,'both');
        if isempty(NonlinearVariables)
            YESNO = nnz(getbasematrix(X,0))==0;
        else
            % No linear terms
            YESNO = isempty(LinearTerms);
            % No constant terms
            YESNO = YESNO && (nnz(getbasematrix(X,0))==0);
            % Largest degree+1
            maxdegree = sum(any(CompressedList(find(ismember(CompressedList,NonlinearVariables)),:),1));
            % All same degree
            YESNO = YESNO && all(all(CompressedList(find(ismember(CompressedList,NonlinearVariables)),2:maxdegree)>0));
        end
        
    case 'factorized'
        YESNO = length(X.leftfactors)>0;
    case 'sos'
        YESNO = (X.typeflag==11);  
    case 'kyp'
        YESNO = (X.typeflag==9);
    case 'gkyp'
        YESNO = (X.typeflag==40);
        
    case 'unitary'
        n = size(X.basis,1);
        YESNO = isequal(X.basis,[sparse(n,1,0) speye(n)]);
        
    otherwise
        error('Wrong input argument.');
end

YESNO = full(YESNO);



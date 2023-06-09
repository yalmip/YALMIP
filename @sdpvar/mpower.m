function y = mpower(x,d)

if (numel(d)>1) || (size(x,1) ~= size(x,2))
   error('Inputs must be a scalar and a square matrix. To compute elementwise POWER, use POWER (.^) instead.');
end

if isa(d,'sdpvar')  
    if numel(x) > 1
        if isa(x,'sdpvar')
            error('x^d not support SDPVARMATRIX^SDPVARSCALAR')
        else
            if isnumeric(x)
                [V,D] = eig(x);
                if ~isreal(D)
                    error('Matrix power x^d requires x to have real eigenvalues');
                end                
                D = real(D);
                y = V*diag(diag(D).^d)*inv(V);
                y.extra.createTime = definecreationtime;
                return                
            else
                error('Object class not support in x^d');
            end
        end
    end
    d.conicinfo = [0 0];
    y = power_internal1(d,x);
    if isa(y,'sdpvar')
        y.extra.createTime = definecreationtime;
        y.extra.opname='';
    end
    return
end
    
% Trivial cases
if d==0
    y = eye(x.dim(1),x.dim(2))^0;
    return
end
if d==1
    y = x;
    return
end
if isnan(d)
	disp('You have NaNs in model (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
	error('NaN power makes no sense.');
end
    

% Check for special case norm(x)^2 and abs(x)^2 which many users try to do
if d==2
    if length(x)==1
        base = getbase(x);
        if isequal(base,[0 1])
            if strcmp(x.extra.opname,'norm')
                model = yalmip('extstruct',getvariables(x));
                z = model.arg{1};
                if (isequal(model.arg{2},2) & min(size(z))==1) | isequal(model.arg{2},'fro')
                    z = reshape(model.arg{1},[],1);
                    y = real(z'*z);
                    y.extra.createTime = definecreationtime;
                    y.extra.opname='';
                    return
                else
                    y = graph_square(x);
                    return
                end
            elseif strcmp(x.extra.opname,'abs')
                model = yalmip('extstruct',getvariables(x));
                z = model.arg{1};
                y = (z.*conj(z));
                y.extra.createTime = definecreationtime;
                y.extra.opname='';
                return                
            end
        end
    end
end


% Fractional and negative powers 
if (ceil(d)-d>0) | (d<0)
    if x.dim(1)>1 | x.dim(2)>1
        error('Only scalars can have negative or non-integer powers');
    else
        base = getbase(x);
        if isequal(base,sparse([0 1])) % Simple unit scalar
            [mt,variabletype,hashM,hash] = yalmip('monomtable');
            var = getvariables(x);
          %  hash = randn(size(mt,2),1);
          %  hashM = mt*hash;
            hashV = (mt(var,:)*d)*hash;
            previous_var = find(abs(hashM - hashV) < 1e-20);
            if isempty(previous_var)
                newmt =  mt(getvariables(x),:)*d;
                mt(end+1,:) = newmt;
                yalmip('setmonomtable',mt,[variabletype newvariabletypegen(newmt)]);
                y = recover(size(mt,1));
            else
                y = recover(previous_var);
            end
        elseif  (size(base,2) == 2) & base(1)==0 & base(2) > 0
            % Something like a*t^-d
            y = base(2)^d*recover(getvariables(x))^d;
        else
            % Bummer, something more complex, add an internal equality constraint 
            y = (yalmip('define','mpower_internal',x))^d;           
        end
    end
    y.extra.createTime = definecreationtime;
    y.extra.opname='';
    return
end

% Integer power of matrix
if x.dim(1)>1 | x.dim(2)>1
    switch d
        case 0
            y = 1;
        case 1
            y = x;
        otherwise
            y = x*mpower(x,d-1);
    end
else %Integer power of scalar
    
    base = x.basis;
    if isequal(base,[0 1]) % Unit scalar can be done fast
        [mt,variabletype,hashes,hash] = yalmip('monomtable');     
        var = x.lmi_variables;
        newmt = mt(var,:)*d;
        newhash = newmt*hash;
        previous_var = findhash(hashes , newhash,length(hashes));
        if isempty(previous_var)
            mt = [mt;newmt];
            yalmip('setmonomtable',mt,[variabletype newvariabletypegen(newmt)]);
            y = x;
            y.lmi_variables = size(mt,1);
        else
            y = x;
            y.lmi_variables = previous_var;
        end
            
    else % General scalar
        switch d
            case 0
                y = 1;
            case 1
                y = x;
            case 2
                y = x*x;
            otherwise
                if even(d)
                    z = mpower(x,d/2);
                    y = z*z;
                else
                    y = x*mpower(x,d-1);
                end
        end
    end
    if isa(y,'sdpvar')
        y.extra.createTime = definecreationtime;
        y.extra.opname='';
    end
end




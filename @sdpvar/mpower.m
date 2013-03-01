function y = mpower(x,d)
%MPOWER (overloaded)

% Author Johan Löfberg
% $Id: mpower.m,v 1.27 2009-10-14 07:28:42 joloef Exp $

%Sanity check
if isa(d,'sdpvar')
    y = power_internal1(d,x);
    return
end

x = flush(x);

if prod(size(d))>1
    error('The power must be scalar.');
end
if x.dim(1)~=x.dim(2) 
    error('Matrix must be square.')
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

% Check for special case norm(x)^2 which many users try to do
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
                    return
                end
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
        elseif  (size(base,2) == 2) & base(1)==0
            % Something like a*t^-d
            y = base(2)^d*recover(getvariables(x))^d;
        else
            % Bummer, something more complex, add an internal equality constraint 
            y = (yalmip('define','mpower_internal',x))^d;           
        end
    end
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
    if 0%isequal(base,[0 1]) % Unit scalar can be done fast
        [mt,variabletype] = yalmip('monomtable');
        %var = getvariables(x);
        var = x.lmi_variables;
        if var > size(mt,2)
            nl_var = find(mt(var,:));
            possible = find(any(mt(:,nl_var),2));
%            possible = 1:size(mt,1);
        else
            possible = find(mt(:,var));
        end
        if length(possible)==1
            % Even faster, we don't need to search, cannot have been
            % definded earlier, since only the linear terms is in monom
            % table                        
            newmt = mt(var,:)*d;
            mt = [mt;newmt];
           % mt(end,size(mt,1))=0;
            yalmip('setmonomtable',mt,[variabletype newvariabletypegen(newmt)]);
            y = x;
            y.lmi_variables = size(mt,1);
        else
            hash = randn(size(mt,2),1);
            mt_hash = mt*hash;
            mt_hash = mt_hash(possible);
            previous_var = findhash(mt_hash , (d*mt(var,:))*hash,size(mt_hash,1));
            if isempty(previous_var)
                newmt = mt(var,:)*d;
                %        mt(end+1,:) = newmt;
                mt = [mt;newmt];
                yalmip('setmonomtable',mt,[variabletype newvariabletypegen(newmt)]);
                y = x;
                y.lmi_variables = size(mt,1);
            else
                previous_var = possible(previous_var);
                y = x;
                y.lmi_variables = previous_var;
            end
        end
    else % General scalar
        switch d
            case 0
                y = 1;
            case 1
                y = x;
            otherwise
                if even(d)
                    z = mpower(x,d/2);
                    y = z*z;
                else
                    y = x*mpower(x,d-1);
                end
        end
    end
end




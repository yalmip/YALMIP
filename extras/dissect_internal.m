function varargout = dissect_internal(X);
% DISSECT_INTERNAL Dissect SDP constraint
%
% G = unblkdiag(F) Converts SDP to several smaller SDPs with more variables
%
% See also UNBLKDIAG

if isa(X,'constraint')
    X = lmi(X);
end

switch class(X)
    case 'sdpvar'
        % Get sparsity pattern
        Z=spy(X);

        % Partition as
        % [A1   0  C1
        %  0   A2  C2
        %  C1' C2' E]

        % Find a dissection
        try
            sep = metismex('NodeBisect',Z);           
        catch
            error('You have to install the MATLAB interface to METIS, http://www.cerfacs.fr/algor/Softs/MESHPART/');
        end

        % Indicies for elements in Ai
        s = setdiff(1:length(Z),sep);

        % re-order Ais to get diagonal blocks
        Z  = Z(s,s);
        AB = X(s,s);
        CD = X(s,sep);
        [v,dummy,r,dummy2]=dmperm(Z);
       
        for i = 1:length(r)-1
            A{i} = AB(v(r(i):r(i+1)-1),v(r(i):r(i+1)-1));
            C{i}= CD(v(r(i):r(i+1)-1),:);
        end
        E = X(sep,sep);
        varargout{1} = A;
        varargout{2} = C;
        varargout{3} = E;

    case 'lmi'
        Fnew=([]);
        % decompose trivial block diagonal stuff
        X = unblkdiag(X);
        for i = 1:length(X)
            if is(X(i),'sdp') & length(sdpvar(X(i)))
                if 0
                [A,B,C,D,E]=dissect(sdpvar(X(i)));
                if ~isempty(B)
                    S=sdpvar(size(E,1));
                    S = S.*inversesparsity(E,D,B);
                    Fnew=Fnew+([A C;C' S]>=0)+([E-S D';D B]>=0);
                else
                    Fnew = Fnew + X(i);
                end
                
                else
                [A,C,E]=dissect(sdpvar(X(i)));
                if length(A)>1
                    allS = 0;
                    for i = 1:length(A)-1
                        S{i}=sdpvar(size(E,1));
                        S{i} = S{i}.*inversesparsity(E,C{i},A{i});
                        allS = allS + S{i};
                        Fnew=Fnew+([A{i} C{i};C{i}' S{i}]>=0);
                    end
                    i = i + 1;
                    S{i}=E-allS;
                    S{i} = S{i}.*inversesparsity(E,C{i},A{i});
                    Fnew=Fnew+([A{i} C{i};C{i}' S{i}]>=0);                                     
                else
                    Fnew = Fnew + X(i);
                end
                                        
                end
                
            else
                Fnew=Fnew+X(i);
            end
        end
        varargout{1} = Fnew;
end


function S = inversesparsity(E,D,B)

if isa(E,'sdpvar')
    E = spy(E).*randn(size(E));
else
    E = E.*randn(size(E));
end
if isa(D,'sdpvar')
    D = spy(D).*randn(size(D));
else
    D = D.*randn(size(D));
end
if isa(B,'sdpvar')
    B = spy(B).*randn(size(B));
else
    B = B.*randn(size(B));
end
S = E-D'*inv(B)*D;
S = S | S';

    

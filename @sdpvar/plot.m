function Y=plot(varargin)
%PLOT (overloaded)

% Fast version for plotting simple PWA objects
if nargin == 1
    X = varargin{1};
    if isa(varargin{1},'sdpvar')
        if length(X) == 1
            if isequal(full(getbase(X)),[zeros(length(X),1) eye(length(X))])
                extstruct = yalmip('extstruct',getvariables(X));
                if ~isempty(extstruct)
                    if isequal(extstruct.fcn,'pwa_yalmip') | isequal(extstruct.fcn,'pwq_yalmip')%#ok
                        % Maybe some reduction has been performed so we
                        % actually can plot it
                        xarg = extstruct.arg{2};                   
                        switch extstruct.arg{3}
                            case ''
                            otherwise
                                Pn = polytope; Bi = {}; Ci = {};
                                index = extstruct.arg{5};
                                for i = 1:length(extstruct.arg{1})
                                    % Pick out row
                                    for j = 1:length(extstruct.arg{1}{i}.Bi)
                                        extstruct.arg{1}{i}.Bi{j} = extstruct.arg{1}{i}.Bi{j}(index,:);
                                        extstruct.arg{1}{i}.Ci{j} = extstruct.arg{1}{i}.Ci{j}(index,:);
                                    end
                                    if isempty(extstruct.arg{1}{i}.Ai{1})
                                        Pn = [Pn extstruct.arg{1}{i}.Pn];
                                        Bi = cat(2, Bi, extstruct.arg{1}{i}.Bi);
                                        Ci = cat(2, Ci, extstruct.arg{1}{i}.Ci);
                                    else
                                        if nnz([extstruct.arg{1}{i}.Ai{:}]) == 0
                                            Pn = [Pn extstruct.arg{1}{i}.Pn];
                                            Bi = cat(2, Bi, extstruct.arg{1}{i}.Bi);
                                            Ci = cat(2, Ci, extstruct.arg{1}{i}.Ci);
                                        else
                                            hold on
                                            [extstruct.arg{1}{i}.Pn,extstruct.arg{1}{i}.Bi,extstruct.arg{1}{i}.Ci,extstruct.arg{1}{i}.Ai] = reduce_basis(extstruct.arg{1}{i}.Pn,extstruct.arg{1}{i}.Bi,extstruct.arg{1}{i}.Ci,extstruct.arg{1}{i}.Ai,xarg);
                                            Y = mpt_plotPWQ(extstruct.arg{1}{i}.Pn, extstruct.arg{1}{i}.Ai,extstruct.arg{1}{i}.Bi,extstruct.arg{1}{i}.Ci);
                                            hold off
                                        end
                                    end
                                end
                                if ~isempty(Bi),
                                    hold on
                                    [Pn,Bi,Ci] = reduce_basis(Pn,Bi,Ci,[],xarg);
                                    if size(Bi{1},2) > 2
                                        error('Cannot plot high-dimensional PWA functions')
                                    end
                                    mpt_plotPWA(Pn, Bi, Ci);
                                    hold off
                                end
                                drawnow
                                return
                        end
                    end
                end
            end
        else
            for j = 1:length(X)
                plot(extsubsref(X,j));
            end
            return
        end
    end
end

% More complex expression. Get epi-graph model
% project to our variables, and extract defining facets
if nargin == 1
    [p,Bi,Ci,Pn,Pfinal] = pwa(varargin{1});%#ok
elseif isa(varargin{2},'lmi')
    [p,Bi,Ci,Pn,Pfinal] = pwa(varargin{1},varargin{2});%#ok
else
    error('Second argument should be a domain defining SET object.');
end
if nargout>0
    Y = mpt_plotPWA(Pn,Bi,Ci);
else
    mpt_plotPWA(Pn,Bi,Ci);
end

function S = extractrow(S,index)
for i = 1:length(S)
    S{i} = S{i}(index,:);
end

function [Pnnew,Binew,Cinew,Ainew] = reduce_basis(Pn,Bi,Ci,Ai,xarg);
%
if ~isequal(getbase(xarg),[zeros(length(xarg),1) eye(length(xarg))])
    base = getbase(xarg);
    c = base(:,1);
    D = base(:,2:end);
    Pnnew = [];
    Ainew = [];
    Binew = [];
    Cinew = [];
    for i = 1:length(Pn)
        [H,K] = double(Pn(i));
        Phere = polytope(H*D,K-H*c);
        if isfulldim(Phere)
            Pnnew = [Pnnew Phere];            
            if ~isempty(Ai)
                Cinew{end+1} = Ci{i} + Bi{i}*c + c'*Ai{i}*c;
                Binew{end+1} = Bi{i}*D + 2*c'*Ai{i}*D;
                Ainew{end+1} = D'*Ai{i}*Ai{i};
            else
                Cinew{end+1} = Ci{i} + Bi{i}*c;
                Binew{end+1} = Bi{i}*D;
            end
        end
    end
else
    Pnnew = Pn;
    Ainew = Ai;
    Binew = Bi;
    Cinew = Ci;
end


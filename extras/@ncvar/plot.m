function Y=plot(varargin)
%PLOT (overloaded)

% Fast version for plotting simple PWA objects
if nargin == 1
    X = varargin{1};
    if isa(varargin{1},'sdpvar')
        if length(X) == 1
            if isequal(full(getbase(X)),[0 1])
                extstruct = yalmip('extstruct',getvariables(X));
                if ~isempty(extstruct)
                    if isequal(extstruct.fcn,'pwa_yalmip') | isequal(extstruct.fcn,'pwq_yalmip')%#ok
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
                                            mpt_plotPWQ(extstruct.arg{1}{i}.Pn, ...
                                                extstruct.arg{1}{i}.Ai, ...
                                                extstruct.arg{1}{i}.Bi, ...
                                                extstruct.arg{1}{i}.Ci, []);
                                            hold off
                                        end
                                    end
                                end
                                if ~isempty(Bi),
                                    mpt_plotPWA(Pn, Bi, Ci);
                                end
                                drawnow
                                return
                        end
                    end
                end
            end
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

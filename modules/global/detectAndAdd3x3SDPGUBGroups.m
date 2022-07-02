function p = detectAndAdd3x3SDPGUBGroups(p)

if any(p.knapsack.type == 3 | p.knapsack.type == 5)
    newAtmost = {};
    top = startofSDPCone(p.K);
    for sdp = 1:length(p.K.s)
        n = p.K.s(sdp);
        for skip = 0:0%n-3
            top_local = top;
            if skip == 0
                max_ = n-2;
            else
                max_ = 1;
            end
            for corner = 1:max_
                n = p.K.s(sdp);
                i = [top_local+0*n top_local+0*n+1+skip top_local+0*n+2+skip,
                     top_local+(1+skip)*n top_local+(1+skip)*n+1+skip top_local+(1+skip)*n+2+skip,
                     top_local+(2+skip)*n top_local+(2+skip)*n+1+skip top_local+(2+skip)*n+2+skip];
                % Extract 3x3
                F_local = p.F_struc(i,:);
                r0 = F_local([1 5 9],2:end);
                if nnz(r0) ==  0
                    % Extract non-diagonal
                    r1 = F_local(2,2:end);
                    r2 = F_local(3,2:end);
                    r3 = F_local(6,2:end);
                    [~,v1,s1] = find(r1);
                    [~,v2,s2] = find(r2);
                    [~,v3,s3] = find(r3);
                    if ~isempty(s1) && ~isempty(s2) && ~isempty(s3) && all(s1(1)==s1) && all(s2(1)==s2) && all(s3(1)==s3)
                        var1 = find_gub_groups(p,v1);
                        var2 = find_gub_groups(p,v2);
                        var3 = find_gub_groups(p,v3);
                        if ~isempty(var1) && ~isempty(var2) && ~isempty(var3)
                            xtest = zeros(length(p.c),1);
                            xtest(v1(1))=1;
                            xtest(v2(1))=1;
                            xtest(v3(1))=1;
                            Xtest = F_local*[1;xtest];
                            Xtest = reshape(Xtest,3,3);
                            if min(eig(Xtest)) < -1e-10
                                allvars = unique([var1 var2 var3]);
                                newAtmost{end+1} = allvars;
                                % Detected a cardinality constraint
                                p.knapsack.a{end+1} = ones(1,length(allvars));
                                p.knapsack.b{end+1} = 2;
                                p.knapsack.variables{end+1} = allvars;
                                p.knapsack.type(end+1) = 2;
                            end
                        end
                    end
                end
                top_local = top_local + n + 1;
            end
        end
        top = top + p.K.s(sdp)^2;
    end
    for i = 1:length(newAtmost)
        row = zeros(1,length(p.c)+1);
        row(1) = 2;
        row(1+newAtmost{i}) = -1;
        p = addInequality(p,row);
%        p.atmost.groups{end+1} = newAtmost{i};
%        p.atmost.bounds(end+1) = 2;
%        p.atmost.variables = unique([p.atmost.variables newAtmost{i}]);
    end
end
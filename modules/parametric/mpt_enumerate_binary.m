function [enums,Matrices] = mpt_enumerate_binary(Matrices)

binary_var_index = Matrices.binary_var_index;
notbinary_var_index = setdiff(1:Matrices.nu,binary_var_index);

% Detect and extract pure binary equalities. These are used
% to detect SOS constraints, and for pruning
nbin = length(binary_var_index);
only_binary = ~any(Matrices.Aeq(:,notbinary_var_index),2);
Aeq_bin = Matrices.Aeq(find(only_binary),binary_var_index);
%Beq_bin = Beq(find(only_binary),binary_var_index);
beq_bin = Matrices.beq(find(only_binary),:);

% Keep the rest
Matrices.Aeq = Matrices.Aeq(find(~only_binary),:);
Matrices.Beq = Matrices.Beq(find(~only_binary),:);
Matrices.beq = Matrices.beq(find(~only_binary),:);

% Dummy pre-solve to avoid problems
if length(beq_bin)>0
    [ii,jj,kk] = unique([Aeq_bin beq_bin],'rows');
    Aeq_bin = Aeq_bin(jj,:);
    beq_bin = beq_bin(jj,:);
end

% Normalize...
neg_b = find(beq_bin < 0);
Aeq_bin(neg_b,:) = -Aeq_bin(neg_b,:);
beq_bin(neg_b,:) = -beq_bin(neg_b,:);

% Detect and extract simple binary LP constraints
only_binary_lp = find(~any([Matrices.G(:,notbinary_var_index) Matrices.E],2));
Abin = Matrices.G(only_binary_lp,binary_var_index);
bbin = Matrices.W(only_binary_lp);

% Remove pure binary constraints
Matrices.G(only_binary_lp,:) = [];
Matrices.W(only_binary_lp,:) = [];
Matrices.E(only_binary_lp,:) = [];

% if ~isempty(Matrices.binary_var_index)
%     fixed_up = find(Matrices.lb(Matrices.binary_var_index) == 1);
%     fixed_down = find(Matrices.ub(Matrices.binary_var_index) == 0);
%     fixed = [fixed_up fixed_down];
%     if ~isempty(fixed)
%         % Ok, some of the binary are fixed.
%         % Do not enumerate over these
%         % Fix them and correct the constraints that we use 
%         % for enumeration and pruning
%         bbin = bbin - Abin(:,fixed)*Matrices.lb(Matrices.binary_var_index(fixed));
%         Abin(:,fixed) = [];
%         beq_bin = beq_bin - Aeq_bin(:,fixed)*Matrices.lb(Matrices.binary_var_index(fixed));
%         Aeq_bin(:,fixed) = [];   
%         Matrices.binary_var_index(fixed) = [];
%         binary_var_index = Matrices.binary_var_index;
%         notbinary_var_index = setdiff(1:Matrices.nu,binary_var_index);
%     end        
% end

% Detect groups with constraints sum(d_i) == d_j
SOS = [];
variables_in_sos = [];
for i = 1:size(Aeq_bin,1)
    if beq_bin(i) == 0
        [ix,jx,sx] = find(Aeq_bin(i,:));
        if all(abs(sx) == 1)
            j = find(sx == -1);
            if length(j) == 1 & length(jx)>2
                % Aha, we have sum binary(i) == binary(j)
                j = jx(j);
                jx = setdiff(jx,j);
                this_sos = sparse(1:length(jx),jx,1,length(jx),length(binary_var_index));
                this_sos(1:end,j)= 1;
                this_sos(end+1,1)=0;
                if isempty(SOS)
                    SOS = this_sos;
                else
%                    new_sos  = [];
                    SOS = kron(ones(size(SOS,1),1),this_sos) | kron(SOS,ones(size(this_sos,1),1));
%                     for k = 1:size(SOS,1)
%                         for r = 1:size(this_sos,1)
%                             new_sos = [new_sos;this_sos(r,:) | SOS(k,:)];
%                         end
%                     end
%                     if ~isequal(new_sos,new_sos2)
%                         error
%                     end
%                    SOS = new_sos;
                end
                variables_in_sos = [variables_in_sos jx j];
            end
        end
    end
end

% Detect groups with constraints sum(d_i) == 1
for i = 1:size(Aeq_bin,1)
    if beq_bin(i) == 1
        [ix,jx,sx] = find(Aeq_bin(i,:));
        if all(sx == 1)
            % Aha, we have sum binary(jx) == 1
            this_sos = sparse(1:length(jx),jx,1,length(jx),length(binary_var_index));
            if isempty(SOS)
                SOS = this_sos;
            else
               % new_sos  = [];
                SOS = kron(ones(size(SOS,1),1),this_sos) | kron(SOS,ones(size(this_sos,1),1));
%                 for k = 1:size(SOS,1)
%                     for r = 1:size(this_sos,1)
%                         new_sos = [new_sos;this_sos(r,:) | SOS(k,:)];
%                     end
%                 end
%                 SOS = new_sos;
            end
            variables_in_sos = [variables_in_sos jx];
        end
    end
end

variables_not_in_sos = setdiff(1:length(binary_var_index),variables_in_sos);
if ~isempty(variables_not_in_sos)
    % we need to add some more binaries
    n_left = length(variables_not_in_sos);
    all_perms = dec2decbin(0:2^n_left-1,n_left);
    [ix,jx,sx] = find(all_perms);jx = variables_not_in_sos(jx);
    this_sos = sparse(ix,jx,sx,size(all_perms,1),length(binary_var_index));
    if isempty(SOS)
        SOS = this_sos;
    else
        new_sos  = [];
        for k = 1:size(SOS,1)
            for r = 1:size(this_sos,1)
                new_sos = [new_sos;this_sos(r,:) | SOS(k,:)];
            end
        end
        SOS = new_sos;
    end
end

enums = unique(SOS,'rows')';

% FIX :FULL, F**N Linux 6.5 repmat bug
feasible = find(~any(Aeq_bin*enums-repmat(full(beq_bin),1,size(enums,2)),1));
enums = enums(:,feasible);

% Remove binary vertices violating these
if ~isempty(Abin)
    remove = find(any(Abin*enums - repmat(bbin,1,size(enums,2)) >0,1));
%    keep = setdiff(1:size(enums,2),remove);
%    mapper = find(ismember(1:size(enums,2),keep));
%    keep(remove)=0;keep2 = 1:size(enums,2);keep2(find(keep)) = find(keep);
%     map = find(keep);
%     remove = find(ismember(j,remove));
%     i(remove) = [];
%     j(remove) = [];
%     k(remove) = [];
%     enums = sparse(i,
%     
%     [i,j,k] = find(enums);
    enums(:,remove)=[];
end

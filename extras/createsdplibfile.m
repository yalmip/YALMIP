function createsdplibfile(F_struc,K, c, filename,integer_variables)
%CREATESDPLIBFILE Internal function to create SDPA data files

fid = fopen(filename,'w');

F_blksz = [];
if any(K.l)
    F_blksz = ones(1,K.l);
else
    F_blksz = [];
end
if any(K.s)
    F_blksz = [F_blksz K.s];
end

try
    fprintf(fid,'%s',['* SDPA format file. Generated by YALMIP ' datestr(now,0)]);fprintf(fid,'\r\n');
    
    n_constraint_matrices = size(F_struc,2)-1;
    n_blocks = size(F_blksz,2);
    
    fprintf(fid,'%i',n_constraint_matrices);fprintf(fid,'\r\n');
    fprintf(fid,'%i',n_blocks);fprintf(fid,'\r\n');
    fprintf(fid,'%i ',F_blksz);fprintf(fid,'\r\n');
    fprintf(fid,'%1.15f ',c);fprintf(fid,'\r\n');
    
    % Another syntax in SDPLIB format
    F_struc(:,1)=-F_struc(:,1);
    
    % Save file (vectorization should be possible but does not seem to improve speed)
    accumsize = [0 cumsum(F_blksz.^2)];
    tol = sqrt(eps);
    for Matrices = 0:n_constraint_matrices
        for blk = 1:n_blocks
            for i = 1:F_blksz(blk)
                for j = i:F_blksz(blk)
                    el = F_struc((i-1)*F_blksz(blk)+j+accumsize(blk),1+Matrices);
                    if abs(el)>tol;
                        fprintf(fid,'%i %i %i %i %1.25f',full([Matrices blk i j el]));fprintf(fid,'\r\n');
                    end
                end
            end
        end
    end
    if ~isempty(integer_variables)
        fprintf(fid,'*INTEGER\r\n');
        for i = 1:length(integer_variables)
            fprintf(fid,'*%i ',integer_variables(i));
        end
    end
    fclose(fid);
catch
    fclose(fid);
end



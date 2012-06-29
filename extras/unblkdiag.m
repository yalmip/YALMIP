function sys = unblkdiag(F)
% UNBLOCKDIAG Extracts diagonal blocks
%
% G = unblkdiag(F) Detects and converts block diagonal SDP.
%
% See also DISSECT

% Author Johan Löfberg
% $Id: unblkdiag.m,v 1.4 2006-01-17 15:49:08 joloef Exp $

switch class(F)
    case 'lmi'
        sys = set([]);
        for i = 1:length(F)
            if is(F(i),'sdp') % SDP
                Z = sdpvar(F(i));
                X = spy(Z);
                [v,dummy,r,dummy2]=dmperm(X);
                if v==dummy & length(r)>2
                    linearblocks = []; % Simple diagonal terms;
                    for blocks = 1:length(r)-1
                        r1 = r(blocks);
                        r2 = r(blocks+1)-1;
                        if r2>r1
                            sys = sys + set(Z(v(r1:r2),v(r1:r2)));
                        else
                            linearblocks = [linearblocks v(r1)];
                        end
                    end
                    if ~isempty(linearblocks)
                        D=diag(Z);
                        sys = sys+set(D(linearblocks));
                    end
                else
                    sys = sys + F(i);
                end
            else
                sys = sys + F(i);
            end
        end
    case 'sdpvar'
        sys = {};
        Z = F;
        X = spy(Z);
        [v,dummy,r,dummy2]=dmperm(X);
        if v==dummy & length(r)>2
            for blocks = 1:length(r)-1
                r1 = r(blocks);
                r2 = r(blocks+1)-1;
                sys{end+1} = Z(v(r1:r2),v(r1:r2));
            end
        else
            sys{1} = Z;
        end
    otherwise
        error('Heh')
end
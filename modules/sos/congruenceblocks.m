function N = congruenceblocks(exponent_m,exponent_p,options,csclasses)
%CONGRUENCEBLOCKS Partitions monomials based on sign symmetry
%
% V = CONGRUENCEBLOCKS(P)
%
% Input
%  V : Vector with SDPVAR objects
%  P : Scalar SDPVAR object
%
% Output
%  V : Cell with SDPVAR objects
%
% Example:
%
% sdpvar x y z
% p = 1+x*y+x^4+y^4+z+z^6;
% v = newtonmonoms(p);
% v = congruenceblocks(v,p);
% sdisplay(v{1}) % Even w.r.t (x,y)
% sdisplay(v{2}) % Odd w.r.t (x,y)
%
% See also NEWTONREDUCE, NEWTONMONOMS, CONSISTENT

% Author Johan Löfberg
% $Id: congruenceblocks.m,v 1.2 2008-11-11 13:29:20 joloef Exp $

sdpvarout = 0;
if isa(exponent_m,'sdpvar')
    z = depends(exponent_p);
    z = recover(unique([depends(exponent_p) depends(exponent_m)]));
    [exponent_p,p_base] = getexponentbase(exponent_p,z);
    [m,m_base] = getexponentbase(exponent_m,z);
    exponent_m = cell(1);exponent_m{1} = m;
    sdpvarout = 1;
end

if nargin < 3
    options.verbose = 0;
    options.sos.congruence = 2;    
end

if nargin < 4
    csclasses = 1;
end

if ~isempty(exponent_m{1}) & options.sos.congruence>0 & ((size(exponent_p,2)<=16)  | options.sos.congruence==1)

    % **********************************************
    % DEFINE CONGRUENCE CLASSES
    % **********************************************
    if options.verbose>0;fprintf('Finding symmetries..............');end;
    n = size(exponent_p,2);
    t = cputime;
    switch options.sos.congruence
        case 1
            Htemp = eye(n);                 % CHEAP VERSION; ONLY CHECK IF IT IS EVEN WRT x_i
        case 2
            Htemp = dec2decbin(1:2^n-1,n)'; % ALL POSSIBLE COMBINATIONS
        otherwise
            error('sos.congruence should be 0, 1 or 2')
    end
    %try
        H = Htemp(:,find(~any(rem(exponent_p*Htemp,2),1))); % Find "even" rows
%     catch
%         i = [];
%         % Loop instead
%         for j = 1:size(Htemp,2)
%             if ~any(rem(exponent_p*Htemp(:,j),2))
%                 i = [i j];
%             end
%         end
%     end
    if isempty(H)
        N = exponent_m;
        if options.verbose>0;disp(['Found no symmetries (' num2str(cputime-t) 'sec)']);end
        return
    end

    t = cputime-t;
    if size(H,2)>=1
        if options.verbose>0
            if size(H,2)>1
                disp(['Found ' num2str(size(H,2)) ' symmetries  (' num2str(t) 'sec)']);
            else
                disp(['Found ' num2str(size(H,2)) ' symmetry  (' num2str(t) 'sec)']);
            end
        end
    else
        if options.verbose>0;disp(['Found no symmetries  (' num2str(t) 'sec)']);end;
    end

    % **********************************************
    % CLASSIFY MONMS ACCORDING TO CONGRUENCE CLASSES
    % **********************************************
    if size(H,2)>=1
        the_text = 'Partitioning using symmetry.....';
    end
    N = cell(0,1);
    for cs = 1:length(csclasses)
        [ur,j,k]=uniquesafe(mod(exponent_m{cs}*H,2),'rows');
        Ntemp = cell(size(ur,1),1);
        temp = [];
        for i = 1:length(k)
            Ntemp{k(i),1} = [Ntemp{k(i)};exponent_m{cs}(i,:)];
            temp = [temp;size(exponent_m{cs}(i,:),1)];
        end

        for i = 1:length(Ntemp)
            N{end+1,1} = Ntemp{i};
        end
    end
    % **********************************************
    % PRINT SOME RESULTS
    % **********************************************
    if size(H,2)>=1
        [uu,ii,oo] = uniquesafe(cellfun('prodofsize',N)/size(N{1},2));
        for i = 1:length(uu)
            n_this = length(find(oo==i));

            the_text = [the_text num2str(uu(i)) 'x' num2str(uu(i)) '(' num2str(n_this) ')' ' '];
        end
    end
    if options.verbose>0;;disp(the_text);end;
else
    N = exponent_m;
end

if sdpvarout
    for i = 1:length(N)
        N{i} = recovermonoms(N{i},z);
    end
end

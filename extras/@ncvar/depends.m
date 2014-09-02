function  LinearVariables = depends(x)
%DEPENDS Returns indicies to variables used in an SDPVAR object
%
% i = depends(x)
%
% Input
%    x : SDPVAR object
% Output
%    i : DOUBLE

[mt,variabletype] = yalmip('monomtable');
ncv = yalmip('nonCommutingTable');

% Simple linear cases
if ~any(variabletype(x.lmi_variables))
    LinearVariables = x.lmi_variables;
else
    LinearVariables = [];
    for i = 1:length(x.lmi_variables)
        v = x.lmi_variables(i);
        if any(mt(v,:))
            LinearVariables = [LinearVariables find(mt(v,:))];
        else
            LinearVariables = [LinearVariables  ncv(v,1+find(ncv(v,2:end)))];
            if ~isnan(ncv(v,1))
                v = ncv(v,1);
                LinearVariables = [LinearVariables find(mt(v,:))];
            end
        end
    end
end
LinearVariables = unique(LinearVariables);
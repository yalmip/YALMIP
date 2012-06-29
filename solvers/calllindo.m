function output = calllindo(interfacedata)

% Author Johan Löfberg
% $Id: calllindo.m,v 1.6 2006-08-18 11:37:13 joloef Exp $

switch interfacedata.solver.tag

    case {'lindo-NLP'}
        output = calllindo_nlp(interfacedata);
    case {'lindo-MIQP'}
        output = calllindo_miqp(interfacedata);
    otherwise
        error;
end

% function output = calllindo_nlp(interfacedata)
% 
% 
% global MY_LICENSE_FILE
% lindo
% 
% % Retrieve needed data
% options = interfacedata.options;
% F_struc = interfacedata.F_struc;
% c       = interfacedata.c;
% K       = interfacedata.K;
% x0      = interfacedata.x0;
% Q       = interfacedata.Q;
% lb      = interfacedata.lb;
% ub      = interfacedata.ub;
% monomtable = interfacedata.monomtable;
% 
% lindo;
% 
% nonlinearindicies = find(interfacedata.variabletype~=0);
% linearindicies = find(interfacedata.variabletype==0);
% nonlinearindicies = union(nonlinearindicies,interfacedata.evalVariables);
% linearindicies    = setdiff(linearindicies,interfacedata.evalVariables);
% interfacedata.nonlinearindicies = nonlinearindicies;
% interfacedata.linearindicies = linearindicies;
% linear = find(interfacedata.variabletype == 0);
% variabletype = interfacedata.variabletype;
% 
% % Init model size
% m  = K.l + K.f;
% n  = length(c);
% csense = [repmat('E',1,K.f) repmat('L',1,K.l)];
% 
% % Specifying variable types...
% vtype = repmat('C',1,length(c(linear)));
% vtype(interfacedata.integer_variables) = 'I';
% 
% oshift = interfacedata.f;
% if m>0
%     A = -F_struc(:,1+linear);
%     b = full(F_struc(:,1));
%     [Nbegcol,Nlencol,Nrowndx,Nobjcnt,Nobjndx,Apatt] = jacSparsity(interfacedata);
%     A = A.*(~Apatt);
% else
%     A = [];
%     b = [];
% end
% 
% [MY_LICENSE_KEY,nErr] = mxlindo('LSloadLicenseString',MY_LICENSE_FILE);
% [iEnv,nErr]=mxlindo('LScreateEnv',MY_LICENSE_KEY);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(-5); return; end;
% [iModel,nErr]=mxlindo('LScreateModel',iEnv);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% constant_data = setup_fmincon_params(interfacedata);
% constant_data.F_struc = F_struc;
% lindo_fun([],[],[],[],[],[],constant_data);
% [nErr] = mxlindo('LSsetFuncalc', iModel, 'lindo_fun',constant_data);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% [nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_PRINTLEVEL, options.verbose);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% 
% % Set NLP solver
% [nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_SOLVER, LS_NMETHOD_MSW_GRG);
% [nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_MAXLOCALSEARCH, 2);
% 
% % Load the LP portion of  model
% [nErr] = mxlindo('LSXloadLPData', iModel, 1, 0, c(linear), b, csense,sparse(A), lb(linear), ub(linear));
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% 
% nErr = mxlindo('LSloadVarType',iModel,vtype);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% 
% % Load the NLP portion of the model
% [nErr] = mxlindo('LSloadNLPData', iModel, Nbegcol, Nlencol,[], Nrowndx, Nobjcnt,Nobjndx,[]);
% if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
% 
% % Optimize model
% solvertime = clock;
% 
% solver = 2;
% solvertime = clock;
% switch solver%interfacedata.solver.tag
%     case 1
%         [solstat,nErr] = mxlindo('LSsolveMIP', iModel);
%         [x,nErr] = mxlindo('LSgetMIPPrimalSolution',iModel);
%     case 2
%         [solstat,nErr] = mxlindo('LSoptimize', iModel,LS_METHOD_FREE);
%         [x,nErr] = mxlindo('LSgetPrimalSolution',iModel);
%     case 3
%         [solStatus,nErr] = mxlindo('LSsolveGOP', iModel);
%         [x,nErr] = mxlindo('LSgetPrimalSolution',iModel);
%     otherwise
% end
% solstat
% if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
% 
% w = zeros(length(c),1);w(linear) =x;
% y = [];
% 
% [nErr]=mxlindo('LSdeleteEnv',iEnv);
% 
% switch solstat
%     case {LS_STATUS_OPTIMAL,LS_STATUS_BASIC_OPTIMAL,7,8}
%         problem = 0;
%     case {LS_STATUS_INFEASIBLE}
%         problem = 1;
%     case {LS_STATUS_UNBOUNDED}
%         problem = 2;
%     otherwise
%         problem = 11;
% end
% infostr = yalmiperror(problem,'LINDO-QP');
% 
% % Save all data sent to solver?
% if options.savesolverinput
%     solverinput.A = A;
%     solverinput.b = b;
%     solverinput.c = c;
%     solverinput.beq = beq;
%     solverinput.options = options.fmincon;
% else
%     solverinput = [];
% end
% 
% % Save all data from the solver?
% if options.savesolveroutput
%     solveroutput.x = x;
%     solveroutput.fmin = fmin;
%     solveroutput.flag = flag;
%     solveroutput.output=output;
%     solveroutput.lambda=lambda;
% else
%     solveroutput = [];
% end
% 
% % Standard interface
% output.Primal      = w;
% output.Dual        = y;
% output.Slack       = [];
% output.problem     = problem;
% output.infostr     = infostr;
% output.solverinput = solverinput;
% output.solveroutput= solveroutput;
% output.solvertime  = solvertime;
% 
% 
% 
% 
% 
% 
% 
% function [Nbegcol,Nlencol,Nrowndx,Nobjcnt,Nobjndx,cJacobian] = jacSparsity(interfacedata)
% 
% linear = find(interfacedata.variabletype == 0);
% oJacobian = zeros(length(linear),1);
% variabletype = interfacedata.variabletype;
% c = interfacedata.c;
% F_struc = interfacedata.F_struc;
% m = size(interfacedata.F_struc,1);
% 
% for i = 1:length(c)
%     if c(i)
%         if variabletype(i)
%             variables = find(interfacedata.monomtable(i,:));
%             oJacobian(variables) = 1;
%         end
%     end
% end
% cJacobian = zeros(m,length(linear));
% for i = 1:size(F_struc,2)-1
%     for j = 1:size(F_struc,1)
%         if F_struc(j,i+1)
%             if variabletype(i)
%                 variables = find(interfacedata.monomtable(i,:));
%                 cJacobian(j,variables) = 1;
%             end
%         end
%     end
% end
% oJacobian = double(oJacobian | any(interfacedata.Q(linear,linear),2));
% 
% Nbegcol = [];
% Nrowndx = [];
% Nlencol = [];
% top = 0;
% for i = 1:size(cJacobian,2)
%     [ii,jj,kk] = find(cJacobian(:,i));
%     if isempty(ii)
%         Nbegcol = [Nbegcol top];
%         Nlencol = [Nlencol 0];
%     else
%         Nbegcol = [Nbegcol top];
%         Nrowndx = [Nrowndx ii(:)'-1];
%         Nlencol = [Nlencol length(ii)];
%         top = top + length(ii);
%     end
% end
% if  isempty(Nrowndx)
%     Nrowndx = [];
% end
% 
% Nobjndx = find(oJacobian) - 1;
% Nobjcnt = length(Nobjndx);
% if  isempty(Nobjndx)
%     Nobjndx = [];
% end
% 
% 
% function output = returnempty(problem)
% output.Primal      = [];
% output.Dual        = [];
% output.Slack       = [];
% output.problem     = problem;
% output.infostr     = yalmiperror(problem,'LINDO');
% output.solverinput = [];
% output.solveroutput= [];
% output.solvertime  = 0;
% 
% 

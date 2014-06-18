function options = sdpsettings(varargin)
%SDPSETTINGS Create/alter solver options structure.
%
%   OPTIONS = SDPSETTINGS with no input arguments returns
%   setting structure with default values
%c
%   OPTIONS = SDPSETTINGS('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
%   OPTIONS = SDPSETTINGS(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%
%
%   SDPSETTINGS PROPERTIES
%
%   GENERAL
%
%    solver             - Specify solver [''|sdpt3|sedumi|sdpa|pensdp|penbmi|csdp|dsdp|maxdet|lmilab|cdd|cplex|xpress|mosek|nag|quadprog|linprog|bnb|bmibnb|kypd|mpt|none ('')]
%    verbose            - Display-level [0|1|2|...(0)] (0-silent, 1-normal, >1-loud)
%    usex0              - Use the current values obtained from double as initial iterate [0|1 (0)]
%    dimacs             - Compute DIMACS error measures [0|1(0)]
%
%    showprogress       - Show progress of YALMIP (suitable for large problems) [0|1 (0)]
%    cachesolvers       - Check for available solvers only first time solvesdp is called [0|1 (0)]
%    warning            - Shows a warning if a problems occurs when solving a problem (infeasibility, numerical problem etc.) [0|1 (1)]
%    beeponproblem      - Beeps when certain warning/error occurs [ integers -2|-1|1|2|3|4|5|6|7|8|9|10|11]
%    saveduals          - Dual variables are saved in YALMIP [0|1 (1)]
%    saveyalmipmodel    - Keep all data sent to solver interface [0|1 (0)]
%    savesolverinput    - Keep all data sent to solver [0|1 (0)]
%    savesolveroutput   - Keep all data returned from solver [0|1 (0)]
%    removeequalities   - Let YALMIP remove equality constraints [-1|0|1 (0)] (-1:with double inequalities, 0:don't, 1:by QR decomposition, 2:basis from constraints)
%    convertconvexquad  - Convert convex quadratic constraints to second order cones [0|1 (1)]
%    radius             - Add radius constraint on all pimal variables ||x||<radius [double >=0 (inf)]
%    shift              - Add small perturbation to (try to) enforce strict feasibility [double >=0 (0)]
%    relax              - Disregard integrality constraint and/or relax nonlinear terms  [0 | 1 (both) 2 (relax integrality) 3 (relax nonlinear terms) (0)]
%    allowmilp          - Allow introduction of binary variables to model nonlinear operators [0 | 1 (0)]
%    expand             - Expand nonlinear operators [0|1 (1)]. Should always be true except in rare debugging cases.
%    plot               - Options when plotting sets
%
%   SUM-OF-SQUARES
%
%    options.sos, see help solvesos
%
%   BRANCH AND BOUND for mixed integer programs
%
%    options.bnb, see help bnb
%
%   BRANCH AND BOUND for polynomial programs
%
%    options.bmibnb, see help bmibnb
%
%   EXTERNAL SOLVERS
%
%    See solver manuals.
%
%   See also SOLVESDP, SET, SDPVAR

% Author Johan Löfberg (ripped from odeset)
% $Id: sdpsettings.m,v 1.80 2010-04-27 14:25:05 joloef Exp $

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help sdpsettings
    return;
end


Names = {'solver'
    'verbose'
    'showprogress'
    'cachesolvers'
    'savesolveroutput'
    'savesolverinput'
    'saveyalmipmodel'
    'saveduals'
    'removeequalities'
    'convertconvexquad'
    'beeponproblem'
    'warning'
    'dimacs'
    'radius'
    'shift'
    'reduce'
    'dualize'
    'relax'
    'usex0'
    'expand'
    'allowmilp'
    'allownonconvex'
    'savedebug'
    'debug'
    %  'plot.colormap'
    %  'plot.gradcolor'
    %'kkt.primalbounds'
    'kkt.dualbounds'
    'kkt.minnormdual'
    'kkt.licqcut'
    'kkt.dualpresolve.lplift'
    'kkt.dualpresolve.passes'
    'plot.edgecolor'
    'plot.wirestyle'
    'plot.wirecolor'
    'plot.linewidth'
    'plot.shade'
    'plot.waitbar'
    'bilevel.algorithm'
    'bilevel.maxiter'
    'bilevel.compslacktol'
    'bilevel.feastol'
    'bilevel.relgaptol'
    'bilevel.rootcuts'
    'bilevel.outersolver'
    'bilevel.innersolver'
    'bilevel.solvefrp'
    'bnb.branchrule'
    'bnb.method'
    'bnb.verbose'
    'bnb.solver'
    'bnb.uppersolver'
    'bnb.presolve'
    'bnb.inttol'
    'bnb.feastol'
    'bnb.gaptol'
    'bnb.maxiter'
    'bnb.round'
    'bnb.weight'
    'bnb.nodefix'
    'bnb.nodetight'
    'bnb.prunetol'
    'bnb.multiple'
    'bnb.profile'
    'bnb.plotbounds'
    'bnb.ineq2eq'
    'bnb.rounding'
    'bmilin.trust'
    'bmilin.alpha'
    'bmilin.beta'
    'bmilin.solver'
    'bmilin.maxiterls'
    'bmilin.maxiter'
    'bmibnb.branchmethod'
    'bmibnb.branchrule'
    'bmibnb.lpreduce'
    'bmibnb.lowrank'
    'bmibnb.diagonalize'
    'bmibnb.lowersolver'
    'bmibnb.uppersolver'
    'bmibnb.lpsolver'
    'bmibnb.target'
    'bmibnb.lowertarget'
    'bmibnb.vartol'
    'bmibnb.relgaptol'
    'bmibnb.absgaptol'
    'bmibnb.pdtol'
    'bmibnb.eqtol'
    'bmibnb.maxiter'
    'bmibnb.maxtime'
    'bmibnb.roottight'
    'bmibnb.numglobal'
    'bmibnb.localstart'
    'bmibnb.strengthscheme'
    'bmibnb.cut.multipliedequality'
    'bmibnb.cut.convexity'
    'bmibnb.cut.evalvariable'
    'bmibnb.cut.bilinear'  
    'bmibnb.cut.monomial'
    'bmibnb.cut.complementarity'
    'cbc.tolint'
    'cbc.maxiter'
    'cbc.maxnodes'
    'cutsdp.solver'
    'cutsdp.maxiter'
    'cutsdp.cutlimit'
    'cutsdp.feastol'
    'cutsdp.variablebound'
    'cutsdp.recoverdual'
    'cutsdp.nodefix'
    'cutsdp.nodetight'
    'mp.algorithm'
    'mp.unbounded'
    'mp.simplify'
    'mp.presolve'
    'robust.lplp'
    'robust.auxreduce'
    'robust.polya'
    'robust.reducedual'
    'robust.reducesemiexplicit'
    'global.branchmethod'
    'global.branchrule'
    'global.lpreduce'
    'global.roottight'
    'global.lowrank'
    'global.lowersolver'
    'global.uppersolver'
    'global.lpsolver'
    'global.target'
    'global.vartol'
    'global.relgaptol'
    'global.absgaptol'
    'global.pdtol'
    'global.eqtol'
    'global.inttol'
    'global.maxiter'
    'global.maxtime'
    
    'gurobi.BarIterLimit'
    'gurobi.Cutoff'
    'gurobi.IterationLimit'
    'gurobi.NodeLimit'
    'gurobi.SolutionLimit'
    'gurobi.TimeLimit'
    'gurobi.BarConvTol'
    'gurobi.BarQCPConvTol'
    'gurobi.FeasibilityTol'
    'gurobi.IntFeasTol'
    'gurobi.MarkowitzTol'
    'gurobi.MIPGap'
    'gurobi.MIPGapAbs'
    'gurobi.OptimalityTol'
    'gurobi.PSDTol'
    'gurobi.InfUnbdInfo'
    'gurobi.NormAdjust'
    'gurobi.ObjScale'
    'gurobi.PerturbValue'
    'gurobi.Quad'
    'gurobi.ScaleFlag'
    'gurobi.Sifting'
    'gurobi.SiftMethod'
    'gurobi.SimplexPricing'
    'gurobi.BarCorrectors'
    'gurobi.BarHomogeneous'
    'gurobi.BarOrder'
    'gurobi.Crossover'
    'gurobi.CrossoverBasis'
    'gurobi.QCPDual'
    'gurobi.BranchDir'
    'gurobi.ConcurrentMIP'
    'gurobi.Heuristics'
    'gurobi.ImproveStartGap'
    'gurobi.ImproveStartNodes'
    'gurobi.ImproveStartTime'
    'gurobi.MinRelNodes'
    'gurobi.MIPFocus'
    'gurobi.MIQCPMethod'
    'gurobi.NodefileDir'
    'gurobi.NodefileStart'
    'gurobi.NodeMethod'
    'gurobi.PumpPasses'
    'gurobi.RINS'
    'gurobi.SolutionNumber'
    'gurobi.SubMIPNodes'
    'gurobi.Symmetry'
    'gurobi.VarBranch'
    'gurobi.ZeroObjNodes'
    'gurobi.TuneOutput'
    'gurobi.TuneResults'
    'gurobi.TuneTimeLimit'
    'gurobi.TuneTrials'
    'gurobi.Cuts'
    'gurobi.CliqueCuts'
    'gurobi.CoverCuts'
    'gurobi.FlowCoverCuts'
    'gurobi.FlowPathCuts'
    'gurobi.GUBCoverCuts'
    'gurobi.ImpliedCuts'
    'gurobi.MIPSepCuts'
    'gurobi.MIRCuts'
    'gurobi.ModKCuts'
    'gurobi.NetworkCuts'
    'gurobi.SubMIPCuts'
    'gurobi.ZeroHalfCuts'
    'gurobi.CutAggPasses'
    'gurobi.CutPasses'
    'gurobi.GomoryPasses'
    'gurobi.AggFill'
    'gurobi.Aggregate'
    'gurobi.DisplayInterval'
    'gurobi.DualReductions'
    'gurobi.FeasRelaxBigM'
    'gurobi.IISMethod'
    'gurobi.LogFile'
    'gurobi.Method'
    'gurobi.NonBlocking'
    'gurobi.NumericFocus'
    'gurobi.PreCrush'
    'gurobi.PreDepRow'
    'gurobi.PreDual'
    'gurobi.PrePasses'
    'gurobi.PreMIQPMethod'
    'gurobi.PrePasses'
    'gurobi.PreQLinearize'
    'gurobi.Presolve'
    'gurobi.PreSparsify'
    'gurobi.ResultFile'
    'gurobi.Threads'
    
    'filtersd.maxiter'
    'filtersd.maxtime'
    'filtersd.maxfeval'        
    'filtersd.rho'
    'filtersd.htol'
    'filtersd.rgtol'
         
    'mpcvx.solver'
    'mpcvx.relgaptol'
    'mpcvx.absgaptol'
    'mpcvx.plot'
    'mpcvx.rays'
    'bpmpd.opts'
    'sp.AbsTol'
    'sp.RelTol'
    'sp.nu'
    'sp.tv'
    'sp.Mfactor'
    'sp.NTiters'
    'maxdet.AbsTol'
    'maxdet.RelTol'
    'maxdet.gam'
    'socp.AbsTol'
    'socp.RelTol'
    'socp.target'
    'socp.max_iter'
    'socp.nu'
    'socp.outmode'
    'sdpt3.vers'
    'sdpt3.gam'
    'sdpt3.predcorr'
    'sdpt3.expon'
    'sdpt3.gaptol'
    'sdpt3.inftol'
    'sdpt3.steptol'
    'sdpt3.maxit'
    'sdpt3.stoplevel'
    'sdpt3.sw2PC_tol'
    'sdpt3.use_corrprim'
    'sdpt3.printyes'
    'sdpt3.scale_data'
    'sdpt3.schurfun'
    'sdpt3.schurfun_parms'
    'sdpt3.randnstate'
    'sdpt3.spdensity'
    'sdpt3.rmdepconstr'
    'sdpt3.CACHE_SIZE'
    'sdpt3.LOOP_LEVEL'
    'sdpt3.cachesize'
    'sdpt3.linsys_options'
    'sdpt3.smallblkdim'
    'sedumi.removeequalities'
    'sedumi.alg'
    'sedumi.beta'
    'sedumi.theta'
    'sedumi.free'
    'sedumi.sdp'
    'sedumi.stepdif'
    'sedumi.w'
    'sedumi.mu'
    'sedumi.eps'
    'sedumi.bigeps'
    'sedumi.maxiter'
    'sedumi.vplot'
    'sedumi.stopat'
    'sedumi.denq'
    'sedumi.denf'
    'sedumi.numtol'
    'sedumi.bignumtol'
    'sedumi.numlvl'
    'sedumi.chol.skip'
    'sedumi.chol.abstol'
    'sedumi.chol.canceltol'
    'sedumi.chol.maxu'
    'sedumi.chol.maxuden'
    'sedumi.chol.numlvl'
    'sedumi.chol.deptol'
    'sedumi.chol.ph1tol'
    'sedumi.cg.qprec'
    'sedumi.cg.restol'
    'sedumi.cg.stagtol'
    'sedumi.cg.maxiter'
    'sedumi.cg.refine'
    
    'sparsecolo.SDPsolver'
    'sparsecolo.domain'
    'sparsecolo.range'
    'sparsecolo.EQorLMI'
    
    'sparsepop.relaxOrder'
    'sparsepop.sparseSW'
    'sparsepop.multiCliquesFactor'
    'sparsepop.scalingSW'
    'sparsepop.boundSW'
    'sparsepop.eqTolerance'
    'sparsepop.perturbation'
    'sparsepop.reduceMomentMatSW'
    'sparsepop.complementaritySW'
    'sparsepop.reduceAMatSW'
    'sparsepop.SDPsolver'
    'sparsepop.SDPsolverSW'
    'sparsepop.SDPsolverEpsilon'
    'sparsepop.SDPsolverOutFile'
    'sparsepop.sdpaDataFile'
    'sparsepop.matFile'
    'sparsepop.POPsolver'
    'sparsepop.detailedInfFile'
    'sparsepop.printFileName'
    'sparsepop.printLevel'
    'sparsepop.errorBdIdx'
    'sparsepop.fValueUbd'
    'sparsepop.symbolicMath'
    'sparsepop.mex'      
    'dsdp.r0'
    'dsdp.zbar'
    'dsdp.penalty'
    'dsdp.boundy'
    'dsdp.gaptol'
    'dsdp.maxit'
    'dsdp.steptol'
    'dsdp.inftol'
    'dsdp.dual_bound'
    'dsdp.dobjmin'
    'dsdp.rho'
    'dsdp.dynamicrho'
    'dsdp.cc'
    'dsdp.bigM'
    'dsdp.mu0'
    'dsdp.reuse'
    'dsdp.lp_barrier'
    'dsdp.objectiveconstant'
    'dsdp.matrixfreesize'
    'dsdp.scale_data'
    'dsdp.maxlanczos'
    'dsdp.max_trust_radius'
    'dsdp.maxtrustradius'
    'dsdp.max_mu_reduction'
    'dsdp.max_infeasible_trust_radius'
    'dsdp.max_infeasible_mu_reduction'
    'sdplr.feastol'
    'sdplr.centol'
    'sdplr.dir'
    'sdplr.penfac'
    'sdplr.reduce'
    'sdplr.limit'
    'sdplr.maxrank'
    'sdplr.soln_factored'
    'sdpa.maxIteration'
    'sdpa.epsilonStar'
    'sdpa.lambdaStar'
    'sdpa.omegaStar'
    'sdpa.lowerBound'
    'sdpa.upperBound'
    'sdpa.betaStar'
    'sdpa.betaBar'
    'sdpa.gammaStar'
    'sdpa.epsilonDash'
    'sdpa.isSymmetric'
    'glpk.lpsolver'
    'glpk.scale'
    'glpk.dual'
    'glpk.price'
    'glpk.relax'
    'glpk.tolbnd'
    'glpk.toldj'
    'glpk.tolpiv'
    'glpk.round'
    'glpk.objll'
    'glpk.objul'
    'glpk.itlim'
    'glpk.tmlim'
    'glpk.branch'
    'glpk.btrack'
    'glpk.tolint'
    'glpk.tolobj'
    'glpk.presol'
    'glpk.save'
    'cplex.presol'
    'cplex.niter'
    'cplex.epgap'
    'cplex.epagap'
    'cplex.relobjdif'
    'cplex.objdif'
    'cplex.tilim'
    'cplex.logfile'
    'cplex.param.double'
    'cplex.param.int'
    'qsopt.dual'
    'qsopt.primalprice'
    'qsopt.dualprice'
    'qsopt.scale'
    'qsopt.maxiter'
    'qsopt.maxtime'
    %'xpress.presol'
    %'xpress.niter'
    'lpsolve.scalemode'
    'cdd.method'
    'clp.solver'
    'clp.maxnumiterations'
    'clp.maxnumseconds'
    'clp.primaltolerance'
    'clp.dualtolerance'
    'clp.primalpivot'
    'clp.dualpivot'
    'pensdp.DEF'
    'pensdp.PBM_MAX_ITER'
    'pensdp.UM_MAX_ITER'
    'pensdp.OUTPUT'
    'pensdp.DENSE'
    'pensdp.LS'
    'pensdp.XOUT'
    'pensdp.UOUT'
    'pensdp.U0'
    'pensdp.MU'
    'pensdp.MU2'
    'pensdp.PBM_EPS'
    'pensdp.P_EPS'
    'pensdp.UMIN'
    'pensdp.ALPHA'
    'pensdp.P0'
    'penbmi.DEF'
    'penbmi.PBM_MAX_ITER'
    'penbmi.UM_MAX_ITER'
    'penbmi.OUTPUT'
    'penbmi.DENSE'
    'penbmi.LS'
    'penbmi.XOUT'
    'penbmi.UOUT'
    'penbmi.NWT_SYS_MODE'
    'penbmi.PREC_TYPE'
    'penbmi.DIMACS'
    'penbmi.TR_MODE'
    'penbmi.U0'
    'penbmi.MU'
    'penbmi.MU2'
    'penbmi.PBM_EPS'
    'penbmi.P_EPS'
    'penbmi.UMIN'
    'penbmi.ALPHA'
    'penbmi.P0'
    'penbmi.PEN_UP'
    'penbmi.ALPHA_UP'
    'penbmi.PRECISION'
    'penbmi.PRECISION_2'
    'penbmi.CG_TOL_DIR'
    
    'pennlp.maxit'
    'pennlp.nwtiters'
    'pennlp.hessianmode'
    'pennlp.autoscale'
    'pennlp.convex'
    'pennlp.eqltymode'
    'pennlp.ignoreinit'
    'pennlp.ncmode'
    'pennlp.nwtstopcrit'
    'pennlp.penalty'
    'pennlp.nwtmode'
    'pennlp.prec'
    'pennlp.cmaxnzs'
    'pennlp.autoini'
    'pennlp.ipenup'
    'pennlp.precision'
    'pennlp.uinit'
    'pennlp.pinit'
    'pennlp.alpha'
    'pennlp.mu'
    'pennlp.dpenup'
    'pennlp.peps'
    'pennlp.umin'
    'pennlp.preckkt'
    'pennlp.cgtolmin'
    'pennlp.cgtolup'
    'pennlp.uinitbox'
    'pennlp.uinitnc'
    
    'lmirank.solver'
    'lmirank.maxiter'
    'lmirank.eps'
    'lmirank.itermod'
    'lmilab.reltol'
    'lmilab.maxiter'
    'lmilab.feasradius'
    'lmilab.L'
    
    'vsdp.solver'
    'vsdp.verifiedupper'
    'vsdp.verifiedlower'
    'vsdp.prove_D_infeasible'
    'vsdp.prove_P_infeasible'
    
    'qpip.mu'
    'qpip.method'
    
    'quadprogbb.max_time'
    'quadprogbb.fathom_tol'
    'quadprogbb.tol'
    'quadprogbb.use_quadprog'
    'quadprogbb.use_single_processor'
    
    'lindo.LS_IPARAM_NLP_SOLVER'
    'lindo.LS_IPARAM_NLP_MAXLOCALSEARCH'
    'lindo.METHOD'   
    'sos.model'
    'sos.newton'
    'sos.congruence'
    'sos.scale'
    'sos.numblkdg'
    'sos.postprocess'
    'sos.sparsetol'
    'sos.impsparse'
    'sos.extlp'
    'sos.csp'
    'sos.savedecomposition'
    'sos.clean'
    'sos.traceobj'
    'sos.reuse'
    'sos.inconsistent'
    'sos.separable'
    'csdp.axtol'
    'csdp.atytol'
    'csdp.objtol'
    'csdp.pinftol'
    'csdp.dinftol'
    'csdp.maxiter'
    'csdp.minstepfrac'
    'csdp.maxstepfrac'
    'csdp.minstepp'
    'csdp.minstepd'
    'csdp.usexzgap'
    'csdp.tweakgap'
    
    'sdpnal.tol'
    'sdpnal.sigma'
    'sdpnal.maxiter'
    'sdpnal.maxitersub'
    'sdpnal.AAtsolve'
    'sdpnal.precond'
    'sdpnal.maxitpsqmr'
    'sdpnal.stagnate_check_psqmr'
    'sdpnal.scale_data'
    'sdpnal.plotyes'
    'sdpnal.proximal'
    
    'specsdp.gatol1'
    'specsdp.gatol2'
    'specsdp.niter=20'
    'specsdp.faiseps'
    'specsdp.compleps'
    'specsdp.stateps'
    'specsdp.penfact'
    'specsdp.penmin'
    'specsdp.pendiv'
    'specsdp.peneps'
    'specsdp.Vinit'
    'specsdp.radius'
    'nag.featol'
    'nag.itmax'
    'nag.bigbnd'
    'nag.orthog'
    'kypd.solver'
    'kypd.lyapunovsolver'
    'kypd.reduce'
    'kypd.transform'
    'kypd.rho'
    'kypd.tol'
    'kypd.lowrank'
    'moment.order'
    'moment.blockdiag'
    'moment.solver'
    'moment.refine'
    'moment.extractrank'
    'moment.rceftol'
    'ipopt.tol'
    'ipopt.mu_strategy'
    'ipopt.hessian_approximation'
    
    'logdetppa.tol'
    'logdetppa.sig'
    'logdetppa.maxiter'
    'logdetppa.maxitersub'
    'logdetppa.precond'
    'logdetppa.maxitpsqmr'
    'logdetppa.stagnate_check_psqmr'
    'logdetppa.scale_data'
    'logdetppa.plotyes'
    'logdetppa.use_proximal'
    'logdetppa.switch_alt_newton_tol'
    
    'powersolver'
    };

% Try to run optimset
[Names,quadprognames,quadprogops]     = trytoset('quadprog',Names);
[Names,linprognames,linprogops]       = trytoset('linprog',Names);
[Names,bintprognames,bintprogops]     = trytoset('bintprog',Names);
[Names,fminconnames,fminconops]       = trytoset('fmincon',Names);
[Names,fminsearchnames,fminsearchops] = trytoset('fminsearch',Names);
[Names,lsqnonnegnames,lsqnonnegops]   = trytoset('lsqnonneg',Names);
[Names,lsqlinnames,lsqlinops]         = trytoset('lsqlin',Names);

obsoletenames ={
    'abstol'
    'reltol'
    'accelerator'
    'maxiters'
    'abstol_init'
    'reltol_init'
    'accelerator_init'
    'maxiters_init'
    'tv'
    'mfactor'
    'sedumi.maxradius'
    'sedumi.uselpsolver'
    'sedumi.removeequalities'
    'silent'
    'sos.solvedual'
    'reduce'
    'penbmi.pbm_eps'};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    % YALMIP options
    options.solver = '';
    options.verbose = 1;
    options.debug = 0;
    options.usex0 = 0;
    options.warning = 1;
    options.cachesolvers = 0;  
    options.showprogress = 0;
    options.saveduals = 1;
    options.removeequalities = 0;
    options.savesolveroutput = 0;
    options.savesolverinput  = 0;
    options.saveyalmipmodel  = 0;
    options.convertconvexquad = 1;
    options.radius = inf;
    options.relax = 0;
    options.dualize = 0;
    options.usex0 = 0;
    options.savedebug = 0;
    options.debug = 0;
    options.expand = 1;
    options.allowmilp = 1;
    options.allownonconvex = 1;
    options.shift = 0;
    options.dimacs = 0;
    options.beeponproblem = [-5 -4 -3 -2 -1];
    % options.plot.colormap = 'hsv';
    % options.plot.gradcolor = 0;
    options.plot.edgecolor = 'k';
    options.plot.wirestyle = '-';
    options.plot.wirecolor = 'k';
    options.plot.linewidth = 0.5;
    options.plot.shade = 1;
    options.plot.waitbar = 1;
    
    %options.kkt.primalbounds = 1;
    options.kkt.dualbounds = 1;
    options.kkt.dualpresolve.passes = 1;
    options.kkt.dualpresolve.lplift = 1;
    options.kkt.minnormdual = 0;
    options.kkt.licqcut = 0;
    
    options.sos.model = 0;
    options.sos.newton = 1;
    options.sos.congruence = 2;
    options.sos.scale = 1;
    options.sos.numblkdg = 0;
    options.sos.postprocess = 0;
    options.sos.csp = 0;
    options.sos.extlp = 1;
    options.sos.impsparse = 0;
    options.sos.sparsetol = 1e-5;
    options.sos.inconsistent = 0;
    options.sos.clean = eps;
    options.sos.savedecomposition = 1;
    options.sos.traceobj = 0;
    options.sos.reuse = 1;
    
    options.moment.order = [];
    options.moment.blockdiag = 0;
    options.moment.solver = '';
    options.moment.refine = 0;
    options.moment.extractrank = 0;
    options.moment.rceftol = -1;
    
    options.robust.lplp = 'enumeration';
    options.robust.auxreduce = 'none';
    options.robust.reducedual = 0;
    options.robust.reducesemiexplicit = 0;
    options.robust.polya = nan;
    
    options.bilevel.algorithm = 'internal';
    options.bilevel.maxiter = 1e4;
    options.bilevel.outersolver = '';
    options.bilevel.innersolver = '';
    options.bilevel.rootcuts = 0;
    options.bilevel.solvefrp = 0;
    options.bilevel.relgaptol = 1e-3;
    options.bilevel.feastol = 1e-6;
    options.bilevel.compslacktol = 1e-8;
    
    
    options.bnb.branchrule = 'max';
    options.bnb.method = 'depthbest';
    options.bnb.verbose = 1;
    options.bnb.solver = '';
    options.bnb.uppersolver = 'rounder';
    options.bnb.presolve = 0;
    options.bnb.inttol = 1e-4;
    options.bnb.feastol = 1e-7;
    options.bnb.gaptol = 1e-4;
    options.bnb.weight = [];
    options.bnb.nodetight = 0;
    options.bnb.nodefix = 0;
    options.bnb.ineq2eq = 0;
    options.bnb.plotbounds = 0;
    options.bnb.rounding = {'ceil','floor','round','shifted round','fix'};
    
    options.bnb.round = 1;
    options.bnb.maxiter = 300;
    
    options.bnb.prunetol = 1e-4;
    options.bnb.multiple = 0;
    options.bnb.profile = 0;
    
    try
        options.bpmpd.opts = bpopt;
    catch
        options.bpmpd.opts =[];
    end
    
    
    % Options for global BMI solver
    options.bmibnb.branchmethod = 'best';
    options.bmibnb.branchrule = 'omega';
    options.bmibnb.cut.multipliedequality = 0;
    options.bmibnb.cut.convexity = 0;
    options.bmibnb.cut.evalvariable = 1;
    options.bmibnb.cut.bilinear = 1;
    options.bmibnb.cut.monomial = 1;
    options.bmibnb.cut.complementarity = 1;
    options.bmibnb.sdpcuts = 0;
    options.bmibnb.lpreduce = 1;
    options.bmibnb.lowrank  = 0;
    options.bmibnb.diagonalize  = 1;
    options.bmibnb.lowersolver = '';
    options.bmibnb.uppersolver = '';
    options.bmibnb.lpsolver = '';
    options.bmibnb.target =  -inf;
    options.bmibnb.lowertarget =  inf;
    options.bmibnb.vartol = 1e-3;
    options.bmibnb.relgaptol = 1e-2;
    options.bmibnb.absgaptol = 1e-2;
    options.bmibnb.pdtol = -1e-6;
    options.bmibnb.eqtol = 1e-6;
    options.bmibnb.maxiter = 100;
    options.bmibnb.maxtime = 3600;
    options.bmibnb.roottight = 1;
    options.bmibnb.numglobal = inf;
    options.bmibnb.localstart = 'relaxed';
    options.bmibnb.presolvescheme = [];
    options.bmibnb.strengthscheme = [1 2 1 3 1 4 1 6 1 5 1 4 1 6 1 4 1];
    
    options.cutsdp.solver = '';
    options.cutsdp.maxiter = 100;
    options.cutsdp.cutlimit = inf;
    options.cutsdp.feastol = -1e-8;
    options.cutsdp.recoverdual = 0;
    options.cutsdp.variablebound = inf;
    options.cutsdp.nodefix = 0;
    options.cutsdp.nodetight = 0;
    
    options.mp.algorithm = 1;
    options.mp.simplify  = 0;
    options.mp.presolve  = 0;
    options.mp.unbounded = 0;
    
    options.global.branchmethod = 'best';
    options.global.branchrule = 'omega';
    options.global.lowersolver = '';
    options.global.uppersolver = '';
    options.global.lpsolver = '';
    options.global.sdpcuts = 0;
    options.global.lpreduce = 1;
    options.global.roottight = 1;
    options.global.lowrank  = 0;
    options.global.target = -inf;
    options.global.vartol = 1e-3;
    options.global.relgaptol = 1e-2;
    options.global.absgaptol = 1e-2;
    options.global.pdtol = -1e-6;
    options.global.eqtol = 1e-6;
    options.global.inttol = 1e-4;
    options.global.round = 1;
    options.global.maxiter = 100;
    options.global.maxtime = 3600;
    
    % Options for approximate multi-parametric
    options.mpcvx.solver = '';
    options.mpcvx.absgaptol = 0.25;
    options.mpcvx.relgaptol = 0.01;
    options.mpcvx.plot = 0;
    options.mpcvx.rays = 'n*20';
    
    options.cbc.tolint = 1e-4;
    options.cbc.maxiter = 10000;
    options.cbc.maxnodes = 100000;

    options.cdd.method = 'criss-cross';
    
    options.clp.solver = 1;
    options.clp.maxnumiterations = 99999999;
    options.clp.maxnumseconds = 3600;
    options.clp.primaltolerance  = 1e-7;
    options.clp.dualtolerance    = 1e-7;
    options.clp.primalpivot = 1;
    options.clp.dualpivot = 1;
    
    options.quadprogbb.max_time = inf;
  
    try
        options.ipopt = ipoptset;
        options.ipopt.hessian_approximation = 'limited-memory';
        
        cNames = recursivefieldnames(options.ipopt);
        for i = 1:length(cNames)
            Names{end+1} = ['ipopt.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);               
    catch
        options.ipopt.mu_strategy = 'adaptive';
        options.ipopt.tol = 1e-7;
        options.ipopt.hessian_approximation = 'limited-memory';    
    end
      
    try   
        options.bonmin = bonminset([],'noIpopt');       
        options.bonmin = rmfield(options.bonmin,'var_lin');
        options.bonmin = rmfield(options.bonmin,'cons_lin');
        
        cNames = recursivefieldnames(options.bonmin);
        for i = 1:length(cNames)
            Names{end+1} = ['bonmin.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.bonmin =[];
    end
    
     try
        options.nomad = nomadset;                       
        cNames = recursivefieldnames(options.nomad);
        for i = 1:length(cNames)
            Names{end+1} = ['nomad.' cNames{i}];          
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.nomad =[];
     end
    
     
     try
         options.clp = clpset;
         cNames = recursivefieldnames(options.clp);
         for i = 1:length(cNames)
             Names{end+1} = ['clp.' cNames{i}];
         end
         [m,n] = size(Names);
         names = lower(Names);
     catch
         
     end
     
     try
         options.cbc = cbcset;
         cNames = recursivefieldnames(options.cbc);
         for i = 1:length(cNames)
             Names{end+1} = ['cbc.' cNames{i}];
         end
         [m,n] = size(Names);
         names = lower(Names);
     catch
         
     end
     
     try
         options.intlinprog =optimoptions('intlinprog');
         evalc('cNames = recursivefieldnames(struct(options.intlinprog));');
         for i = 1:length(cNames)
             Names{end+1} = ['intlinprog.' cNames{i}];
         end
         [m,n] = size(Names);
         names = lower(Names);
     catch
         options.intlinprog = [];
     end
     
     try
         options.ooqp = ooqpset;
         cNames = recursivefieldnames(options.ooqp);
         for i = 1:length(cNames)
             Names{end+1} = ['ooqp.' cNames{i}];
         end
         [m,n] = size(Names);
         names = lower(Names);
     catch
         options.ooqp = [];
     end
    
    try
        options.xpress = xprsoptimset;
        cNames = recursivefieldnames(options.xpress);
        for i = 1:length(cNames)
            Names{end+1} = ['xpress.' cNames{i}];
            options.xpress = setfield(options.xpress,cNames{i},[]);
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.xpress =[];
    end
            
    try
        options.qpoases = qpOASES_options;
        cNames = recursivefieldnames(options.qpoases);
        for i = 1:length(cNames)
            Names{end+1} = ['xpress.' cNames{i}];
            options.qpoases = setfield(options.qpoases,cNames{i},[]);
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.qpoases =[];
    end
        
    try
        options.cplex = cplexoptimset('cplex');
        options.cplex.output.clonelog = 0;
        remove = zeros(length(Names),1);
        for i = 1:length(Names)
            if strfind(Names{i},'cplex')
                remove(i)=1;
            end
        end
        Names = Names(find(~remove));
        cNames = recursivefieldnames(options.cplex);
        for i = 1:length(cNames)
            Names{end+1} = ['cplex.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.cplex.presol = 1;
        options.cplex.niter = 1;
        options.cplex.epgap = 1e-4;
        options.cplex.epagap = 1e-6;
        options.cplex.relobjdif = 0.0;
        options.cplex.objdif = 0.0;
        options.cplex.tilim = 1e75;
        options.cplex.logfile = 0;
        options.cplex.param.double = [];
        options.cplex.param.int = [];
    end
    
    try
        options.baron = baronset;
        cNames = recursivefieldnames(options.baron);
        for i = 1:length(cNames)
            Names{end+1} = ['baron.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.baron = [];
    end
    
    try
        options.knitro = optimset;
        options.knitro.optionsfile = '';
        cNames = recursivefieldnames(options.knitro);
        for i = 1:length(cNames)
            Names{end+1} = ['knitro.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
       options.knitro.optionsfile = '';
    end
    
    try
        % OPTI Toolbox interface
        options.csdp = csdpset();
        
        remove = zeros(length(Names),1);
        for i = 1:length(Names)
            if strfind(Names{i},'csdp')
                remove(i)=1;
            end
        end
        Names = Names(find(~remove));
        cNames = recursivefieldnames(options.csdp);
        for i = 1:length(cNames)
            Names{end+1} = ['csdp.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
        
    catch
        options.csdp.axtol  = 1e-8;
        options.csdp.atytol = 1e-8;
        options.csdp.objtol = 1e-8;
        options.csdp.pinftol = 1e8;
        options.csdp.dinftol = 1e8;
        options.csdp.maxiter = 100;
        options.csdp.minstepfrac = 0.90;
        options.csdp.maxstepfrac = 0.97;
        options.csdp.minstepp = 1e-8;
        options.csdp.minstepd = 1e-8;
        options.csdp.usexzgap = 1;
        options.csdp.tweakgap = 0;
    end
        
    try
        options.scip = optiset;
        cNames = recursivefieldnames(options.scip);
        for i = 1:length(cNames)
            Names{end+1} = ['scip.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.scip = [];
    end
    
    try
        % OPTI Toolbox interface
        options.dsdp = dsdpset();
        
        remove = zeros(length(Names),1);
        for i = 1:length(Names)
            if strfind(Names{i},'dsdp')
                remove(i)=1;
            end
        end
        Names = Names(find(~remove));
        cNames = recursivefieldnames(options.dsdp);
        for i = 1:length(cNames)
            Names{end+1} = ['dsdp.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
        
    catch
        % Options for DSDP 5.6 classical interface
        options.dsdp.r0 = -1;
        options.dsdp.zbar = 0;
        options.dsdp.penalty  = 1e8;
        options.dsdp.boundy  = 1e6;
        options.dsdp.gaptol = 1e-7;
        options.dsdp.maxit  = 500;
        options.dsdp.steptol=5e-2;
        options.dsdp.inftol=1e-8;
        options.dsdp.dual_bound = 1e20;
        options.dsdp.rho = 3;
        options.dsdp.dynamicrho = 1;
        options.dsdp.bigM = 0;
        options.dsdp.mu0 = -1;
        options.dsdp.reuse = 4;
        options.dsdp.lp_barrier = 1;
    end
    
    try        
        evalc('[r,res]=mosekopt(''param'');');
        options.mosek = res.param;
        remove = zeros(length(Names),1);
        for i = 1:length(Names)
            if strfind(Names{i},'mosek')
                remove(i)=1;
            end
        end
        Names = Names(find(~remove));
        cNames = recursivefieldnames(options.mosek);
        for i = 1:length(cNames)
            Names{end+1} = ['mosek.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.mosek.param = [];
    end
    
    try
        options.penlab = penlab.defopts(1);
        cNames = recursivefieldnames(options.penlab);
        for i = 1:length(cNames)
            Names{end+1} = ['penlab.' cNames{i}];
        end
        [m,n] = size(Names);
        names = lower(Names);
    catch
        options.penlab = [];
    end
    
    
%     options.csdp.axtol  = 1e-8;
%     options.csdp.atytol = 1e-8;
%     options.csdp.objtol = 1e-8;
%     options.csdp.pinftol = 1e8;
%     options.csdp.dinftol = 1e8;
%     options.csdp.maxiter = 100;
%     options.csdp.minstepfrac = 0.90;
%     options.csdp.maxstepfrac = 0.97;
%     options.csdp.minstepp = 1e-8;
%     options.csdp.minstepd = 1e-8;
%     options.csdp.usexzgap = 1;
%     options.csdp.tweakgap = 0;    
  %  options.mosek.param = [];
    
%     % Options for DSDP 5.6
%     options.dsdp.r0 = -1;
%     options.dsdp.zbar = 0;
%     options.dsdp.penalty  = 1e8;
%     options.dsdp.boundy  = 1e6;
%     options.dsdp.gaptol = 1e-7;
%     options.dsdp.maxit  = 500;
%     options.dsdp.steptol=5e-2;
%     options.dsdp.inftol=1e-8;
%     options.dsdp.dual_bound = 1e20;
%     options.dsdp.rho = 3;
%     options.dsdp.dynamicrho = 1;
%     options.dsdp.bigM = 0;
%     options.dsdp.mu0 = -1;
%     options.dsdp.reuse = 4;
%     options.dsdp.lp_barrier = 1;
    
    % Older versions
    options.dsdp.objectiveconstant = 0;
    options.dsdp.matrixfreesize = 3000;
    options.dsdp.scale_data = 1;
    options.dsdp.max_trust_radius = -1;
    options.dsdp.maxtrustradius = 1e10;
    options.dsdp.max_infeasible_trust_radius = -1;
    options.dsdp.max_infeasible_mu_reduction = 2;
    options.dsdp.max_mu_reduction = 1e8;
    options.dsdp.maxlanczos = 20;
    
    
    options.filtersd.maxiter = 1500;
    options.filtersd.maxtime = 1000;
    options.filtersd.maxfeval = 10000;         
    
    % Options for GLPK
    options.glpk.lpsolver = 1;
    options.glpk.scale = 1;
    options.glpk.dual = 0;
    options.glpk.price = 1;
    options.glpk.relax = 0.07;
    options.glpk.tolbnd = 1e-7;
    options.glpk.toldj = 1e-7;
    options.glpk.tolpiv = 1e-9;
    options.glpk.round = 0;
    options.glpk.objll = -1e12;
    options.glpk.objul = 1e12;
    options.glpk.itlim = 1e4;
    options.glpk.tmlim = -1;
    options.glpk.branch = 2;
    options.glpk.btrack = 2;
    options.glpk.tolint = 1e-6;
    options.glpk.tolobj = 1e-7;
    options.glpk.presol = 1;
    options.glpk.save = 0;
            
    options.gurobi.BarIterLimit = inf;
    options.gurobi.Cutoff = inf;
    options.gurobi.IterationLimit = inf;
    options.gurobi.NodeLimit = inf;
    options.gurobi.SolutionLimit = inf;
    options.gurobi.TimeLimit = inf;
    
    options.gurobi.BarConvTol = 1e-8;
    options.gurobi.BarQCPConvTol = 1e-6;
    options.gurobi.FeasibilityTol = 1e-6;
    options.gurobi.IntFeasTol = 1e-6;
    options.gurobi.MarkowitzTol = 0.0078125;
    options.gurobi.MIPGap = 1e-4;
    options.gurobi.MIPGapAbs = 1e-10;
    options.gurobi.OptimalityTol = 1e-6;
    options.gurobi.PSDTol = 1e-6;
    options.gurobi.InfUnbdInfo = 0;
    options.gurobi.NormAdjust = -1;
    options.gurobi.ObjScale = 0;
    options.gurobi.PerturbValue = 0.0002;
    options.gurobi.Quad = -1;
    options.gurobi.ScaleFlag = 1;
    options.gurobi.Sifting = -1;
    options.gurobi.SiftMethod = -1;
    options.gurobi.SimplexPricing = -1;
    options.gurobi.BarCorrectors = -1;
    options.gurobi.BarHomogeneous = -1;
    options.gurobi.BarOrder = -1;
    options.gurobi.Crossover = -1;
    options.gurobi.CrossoverBasis = 0;
    options.gurobi.QCPDual = 0;
    options.gurobi.BranchDir = 0;
    options.gurobi.ConcurrentMIP = 1;
    options.gurobi.Heuristics = 0.05;
    options.gurobi.ImproveStartGap = 0;
    options.gurobi.ImproveStartNodes = inf;
    options.gurobi.ImproveStartTime = inf;
    options.gurobi.MinRelNodes = 0;
    options.gurobi.MIPFocus = 0;
    options.gurobi.MIQCPMethod = -1;
    options.gurobi.NodefileDir = '.';
    options.gurobi.NodefileStart = inf;
    options.gurobi.NodeMethod = 1;
    options.gurobi.PumpPasses = 0;
    options.gurobi.RINS = -1;
    options.gurobi.SolutionNumber = 0;
    options.gurobi.SubMIPNodes = 500;
    options.gurobi.Symmetry = -1;
    options.gurobi.VarBranch = -1;
    options.gurobi.ZeroObjNodes = 0;
    options.gurobi.TuneOutput = 2;
    options.gurobi.TuneResults = -1;
    options.gurobi.TuneTimeLimit = -1;
    options.gurobi.TuneTrials = 2;
    options.gurobi.Cuts = -1;
    options.gurobi.CliqueCuts = -1;
    options.gurobi.CoverCuts = -1;
    options.gurobi.FlowCoverCuts = -1;
    options.gurobi.FlowPathCuts = -1;
    options.gurobi.GUBCoverCuts = -1;
    options.gurobi.ImpliedCuts = -1;
    options.gurobi.MIPSepCuts = -1;
    options.gurobi.MIRCuts = -1;
    options.gurobi.ModKCuts = -1;
    options.gurobi.NetworkCuts = -1;
    options.gurobi.SubMIPCuts = -1;
    options.gurobi.ZeroHalfCuts = -1;
    options.gurobi.CutAggPasses = -1;
    options.gurobi.CutPasses = -1;
    options.gurobi.GomoryPasses = -1;
    options.gurobi.AggFill = 10;
    options.gurobi.Aggregate = 1;    
    options.gurobi.DisplayInterval = 5;
    options.gurobi.DualReductions = 1;
    options.gurobi.FeasRelaxBigM = 1e6;
    options.gurobi.IISMethod = -1;
    options.gurobi.LogFile = '';
    options.gurobi.Method = -1;
    options.gurobi.NonBlocking = 0;
    options.gurobi.NumericFocus = 0;
    options.gurobi.PreCrush = 0;
    options.gurobi.PreDepRow = -1;
    options.gurobi.PreDual = -1;
    options.gurobi.PreMIQPMethod = -1;
    options.gurobi.PreQLinearize = -1;
    options.gurobi.PrePasses = -1;
    options.gurobi.Presolve = -1;
    options.gurobi.PreSparsify = 0;
    options.gurobi.ResultFile = '';
    options.gurobi.Threads = 0;
      
    options.kypd.solver = '';
    options.kypd.lyapunovsolver = 'schur';
    options.kypd.reduce = 0;
    options.kypd.transform = 0;
    options.kypd.rho = 1;
    options.kypd.tol = 1e-8;
    options.kypd.lowrank = 0;
    
    options.lmilab.reltol = 1e-3;
    options.lmilab.maxiter = 100;
    options.lmilab.feasradius = 1e9;
    options.lmilab.L = 10;
    
    options.lmirank.solver = '';
    options.lmirank.maxiter = 100;
    options.lmirank.maxiter = 1000;
    options.lmirank.eps = 1e-9;
    options.lmirank.itermod = 1;
    
    
    options.logdetppa.tol = 1e-6;
    options.logdetppa.sig = 10;
    options.logdetppa.maxiter    = 100;
    options.logdetppa.maxitersub = 30;
    options.logdetppa.precond    = 1;
    options.logdetppa.maxitpsqmr = 100;
    options.logdetppa.stagnate_check_psqmr = 0;
    options.logdetppa.scale_data = 2;
    options.logdetppa.plotyes    = 0;
    options.logdetppa.use_proximal = 1;
    options.logdetppa.switch_alt_newton_tol = 1e-2;
    
    options.vsdp.solver = '';
    options.vsdp.verifiedupper = 0;
    options.vsdp.verifiedlower = 1;
    options.vsdp.prove_D_infeasible = 0;
    options.vsdp.prove_P_infeasible = 0;
    
    options.lindo.LS_IPARAM_NLP_SOLVER = 'LS_NMETHOD_MSW_GRG';
    options.lindo.LS_IPARAM_NLP_MAXLOCALSEARCH = 1;
    options.lindo.LS_METHOD = 'LS_METHOD_FREE';
    
    options.lpsolve.scalemode = 0;
    
    % Options for MAXDET solver
    options.maxdet.AbsTol = 1e-6;
    options.maxdet.RelTol = 1e-6;
    options.maxdet.gam    = 25;
    options.maxdet.NTiters= 50;
    
    options.nag.featol = sqrt(eps);
    options.nag.itmax = 500;
    options.nag.bigbnd = 1e10;
    options.nag.orthog = 0;
    
    options.penbmi.DEF = 1;
    options.penbmi.PBM_MAX_ITER = 50;
    options.penbmi.UM_MAX_ITER = 100;
    options.penbmi.OUTPUT = 1;
    options.penbmi.DENSE = 1;          %!0
    options.penbmi.LS = 0;
    options.penbmi.XOUT = 0;
    options.penbmi.UOUT = 0;
    options.penbmi.NWT_SYS_MODE = 0;
    options.penbmi.PREC_TYPE = 0;
    options.penbmi.DIMACS = 0;
    options.penbmi.TR_MODE = 0;
    options.penbmi.U0 = 1;
    options.penbmi.MU = 0.7;
    options.penbmi.MU2 = 0.5;          %!0.1
    % options.penbmi.PBM_EPS = 1e-6;     %!1e-7
    options.penbmi.PRECISION = 1e-6;     %!1e-7
    options.penbmi.P_EPS = 1e-4;       %!1e-6
    options.penbmi.UMIN = 1e-14;
    options.penbmi.ALPHA = 1e-2;
    options.penbmi.P0 = 0.1;           %!0.01
    options.penbmi.PEN_UP = 0.5;       %!0
    options.penbmi.ALPHA_UP = 1.0;
    options.penbmi.PRECISION_2 = 1e-6; %!1e-7
    options.penbmi.CG_TOL_DIR = 5e-2;
    
    options.pennlp.maxit = 100;
    options.pennlp.nwtiters = 100;
    options.pennlp.hessianmode = 0;
    options.pennlp.autoscale = 1;
    options.pennlp.convex = 0;
    options.pennlp.eqltymode = 1;
    options.pennlp.ignoreinit = 0;
    options.pennlp.ncmode = 0;
    options.pennlp.nwtstopcrit = 2;
    options.pennlp.penalty = 0;
    options.pennlp.nwtmode = 0;
    options.pennlp.prec = 0;
    options.pennlp.cmaxnzs =-1;
    options.pennlp.autoini = 1;
    options.pennlp.ipenup = 1;
    
    options.pennlp.precision = 1e-7;
    options.pennlp.uinit = 1;
    options.pennlp.pinit = 1;
    options.pennlp.alpha = 0.01;
    options.pennlp.mu = 0.5;
    options.pennlp.dpenup = 0.1;
    options.pennlp.peps = 1e-8;
    options.pennlp.umin = 1e-12;
    options.pennlp.preckkt = 1e-1;
    options.pennlp.cgtolmin = 5e-2;
    options.pennlp.cgtolup = 1;
    options.pennlp.uinitbox = 1;
    options.pennlp.uinitnc = 1;
    
    
    options.pensdp.DEF = 1;
    options.pensdp.PBM_MAX_ITER = 50;
    options.pensdp.UM_MAX_ITER = 100;
    options.pensdp.OUTPUT = 1;
    options.pensdp.DENSE = 0;
    options.pensdp.LS = 0;
    options.pensdp.XOUT = 0;
    options.pensdp.UOUT = 0;
    options.pensdp.U0 = 1;
    options.pensdp.MU = 0.7;
    options.pensdp.MU2 = 0.1;
    options.pensdp.PBM_EPS = 1e-7;
    options.pensdp.P_EPS = 1e-6;
    options.pensdp.UMIN = 1e-14;
    options.pensdp.ALPHA = 1e-2;
    options.pensdp.P0 = 0.9;
    
    options.powersolver = [];
    
    
    options.quadprogbb.max_time = 86400;
    options.quadprogbb.fathom_tol = 1e-6;
    options.quadprogbb.tol = 1e-8;
    options.quadprogbb.use_quadprog = 1;
    options.quadprogbb.use_single_processor = 0;
    
    
    
    % Options for SDPA
    options.sdpa.maxIteration = 100 ;
    options.sdpa.epsilonStar = 1.0E-7;
    options.sdpa.lambdaStar  = 1.0E2  ;
    options.sdpa.omegaStar  = 2.0 ;
    options.sdpa.lowerBound  = -1.0E5  ;
    options.sdpa.upperBound  = 1.0E5  ;
    options.sdpa.betaStar  = 0.1  ;
    options.sdpa.betaBar  = 0.2 ;
    options.sdpa.gammaStar  = 0.9 ;
    options.sdpa.epsilonDash  = 1.0E-7 ;
    options.sdpa.isSymmetric = 0 ;
    
    options.sdplr.feastol = 1e-5;
    options.sdplr.centol = 1e-1;
    options.sdplr.dir = 1;
    options.sdplr.penfac = 2;
    options.sdplr.reduce = 0;
    options.sdplr.limit = 3600;
    options.sdplr.soln_factored = 0;
    options.sdplr.maxrank = 0;
    
    % Options for SDPT3 (2.3 -> 3.01)
    options.sdpt3.vers     = 1;
    options.sdpt3.gam      = 0;
    options.sdpt3.predcorr = 1;
    options.sdpt3.expon    = 1;
    options.sdpt3.gaptol   = 1e-7;
    options.sdpt3.inftol   = 1e-7;
    options.sdpt3.steptol  = 1e-6;
    options.sdpt3.maxit    = 50;
    options.sdpt3.stoplevel= 1;
    options.sdpt3.sw2PC_tol  = inf;
    options.sdpt3.use_corrprim  = 0;
    options.sdpt3.printyes   = 1;
    options.sdpt3.scale_data = 0;
    options.sdpt3.schurfun   = [];
    options.sdpt3.schurfun_parms = [];
    options.sdpt3.randnstate = 0;
    options.sdpt3.spdensity   = 0.5;
    options.sdpt3.rmdepconstr = 0;
    options.sdpt3.CACHE_SIZE = 256;
    options.sdpt3.LOOP_LEVEL = 8;
    options.sdpt3.cachesize = 256;
    options.sdpt3.linsys_options = 'raugmatsys';
    options.sdpt3.smallblkdim = 30;
    
    options.sparsecolo.SDPsolver = '';
    options.sparsecolo.domain = 2;
    options.sparsecolo.range = 1;
    options.sparsecolo.EQorLMI = 2;
    
    options.sparsepop.relaxOrder = 1;
    options.sparsepop.sparseSW = 1;
    options.sparsepop.multiCliquesFactor = 1;
    options.sparsepop.scalingSW = 1;
    options.sparsepop.boundSW = 2;
    options.sparsepop.eqTolerance = 0;
    options.sparsepop.perturbation = 0;
    options.sparsepop.reduceMomentMatSW = 1;
    options.sparsepop.complementaritySW = 0;
    options.sparsepop.reduceAMatSW = 1;
    options.sparsepop.SDPsolver = 'sedumi';
    options.sparsepop.SDPsolverSW = 1;
    options.sparsepop.SDPsolverEpsilon = 1.0000e-009;
    options.sparsepop.SDPsolverOutFile = 0;
    options.sparsepop.sdpaDataFile = '';
    options.sparsepop.matFile = '';
    options.sparsepop.POPsolver = '';
    options.sparsepop.detailedInfFile = '';
    options.sparsepop.printFileName = 1;
    %  options.sparsepop.printLevel = [2 0];
    options.sparsepop.errorBdIdx = '';
    options.sparsepop.fValueUbd = '';
    options.sparsepop.symbolicMath = 1;
    options.sparsepop.mex = 0;
    
    options.sdpnal.tol = 1e-6;
    options.sdpnal.sigma = 10;
    options.sdpnal.maxiter = 100;
    options.sdpnal.maxitersub = 20;
    options.sdpnal.AAtsolve = 2;
    options.sdpnal.precond = 1;
    options.sdpnal.maxitpsqmr = 100;
    options.sdpnal.stagnate_check_psqmr = 0;
    options.sdpnal.scale_data = 2;
    options.sdpnal.plotyes = 0;
    options.sdpnal.proximal = 1;
    
    % Options for SeDuMi (1.03 -> 1.05)
    options.sedumi.alg    = 2;
    options.sedumi.beta   = 0.5;
    options.sedumi.theta  = 0.25;
    options.sedumi.free   = 1;
    options.sedumi.sdp    = 0;
    options.sedumi.stepdif= 0;
    options.sedumi.w      = [1 1];
    options.sedumi.mu     = 1.0;
    options.sedumi.eps    = 1e-9;
    options.sedumi.bigeps = 1e-3;
    options.sedumi.maxiter= 150;
    options.sedumi.vplot  = 0;
    options.sedumi.stopat     = -1;
    options.sedumi.denq   = 0.75;
    options.sedumi.denf   = 10;
    options.sedumi.numtol = 5e-7;
    options.sedumi.bignumtol = 0.9;
    options.sedumi.numlvlv = 0;
    options.sedumi.chol.skip = 1;
    options.sedumi.chol.canceltol = 1e-12;
    options.sedumi.chol.maxu   = 5e5;
    options.sedumi.chol.abstol = 1e-20;
    options.sedumi.chol.maxuden= 5e2;
    options.sedumi.cg.maxiter = 25;
    options.sedumi.cg.restol  = 5e-3;
    options.sedumi.cg.refine  = 1;
    options.sedumi.cg.stagtol = 5e-14;
    options.sedumi.cg.qprec   = 0;
    options.sedumi.maxradius = inf;
    
    
    options.qsopt.dual = 0;
    options.qsopt.primalprice = 1;
    options.qsopt.dualprice = 6;
    options.qsopt.scale = 1;
    options.qsopt.maxiter = 300000;
    options.qsopt.maxtime = 10000.0;
    
    options.qpip.mu = 0.0;
    options.qpip.method = 1;
    
    %  options.xpress.presol = 1;
    %  options.xpress.niter  = 1;
    
    try
        % Maybe we already created these above
        if exist('fminconops','var')
            options.quadprog = quadprogops;
            options.linprog  = linprogops;
            options.bintprog  = bintprogops;
            options.fmincon  = fminconops;
            options.fminsearch  = fminsearchops;
            options.lsqnonneg  = lsqnonnegops;
            options.lsqlin  = lsqlinops;
        else
            options.quadprog = optimset('quadprog');
            options.linprog  = optimset('linprog');
            options.bintprog  = optimset('linprog');
            options.fmincon  = optimset('fmincon');
            options.fminsearch  = optimset('fminsearch');
            options.lsqnonneg  = optimset('lsqnonneg');
            options.lsqlin  = optimset('lsqlin');
        end
    catch
        try
            % Ok, what about tomlab then?
            options.quadprog = optimset;
            options.linprog  = optimset;
            options.bintprog = optimset;
            options.fmincon  = optimset;
            options.fminsearch  = optimset;
        catch
            options.quadprog = [];
            options.linprog  = [];
            options.bintprog = [];
            options.fmincon  = [];
            options.fminsearch  = [];
        end
    end
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability...
            
            if strcmp(lowArg,'sedumi.removeequalities')
                lowArg = 'removeequalities';
                arg = lowArg;
                
            elseif strcmp(lowArg,'silent')
                lowArg = 'verbose';
                if varargin{i+1}==1
                    varargin{i+1} = 0;
                else
                    varargin{i+1} = 1;
                end
                warning('The field ''silent'' is obsolete. Use ''verbose'' instead')
            elseif strcmp(lowArg,'sedumi.maxradius')
                lowArg = 'radius';
                arg = lowArg;
                warning('The field ''sedumi.maxradius'' is obsolete. Use ''radius'' instead')
            elseif strcmp(lowArg,'penbmi.pbm_eps')
                lowArg = 'penbmi.precision';
                arg = lowArg;
                warning('The field ''penbmi.pbm_eps'' is obsolete. Use ''penbmi.precision'' instead')
            elseif strcmp(lowArg,'sos.solvedual')
                lowArg = 'sos.model';
                arg = lowArg;
            else
                error(sprintf('The property name ''%s'' is obsolete. Sorry...', arg));
            end
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if (length(k) == 1)
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end


function [Names,solvernames,solverops] = trytoset(solver,Names)

try
    try
        try
            evalc(['solverops = ' solver '(''defaults'');']);
        catch
            solverops = optimset(solver);
        end
    catch
        solverops = optimset;
    end
    solvernames = fieldnames(solverops);
    
    
    if isequal(solver, 'quadprog') && isfield(solverops, 'Algorithm') && ~isempty(solverops.Algorithm)
        solverops.Algorithm = 'active-set';
    end

    if  any(strcmp(solvernames,'LargeScale'))
        if isequal(solver, 'quadprog')
            solverops.LargeScale = 'off';
        end
    else
        solvernames{end+1} = 'LargeScale';
        solverops.LargeScale = 'off';
    end
    
    for i = 1:length(solvernames)
        Names = [Names;{[solver '.' solvernames{i}]}];
    end
    
catch
    solvernames = {'LargeScale'};
    solverops = struct('LargeScale','off');
end


function cNames = recursivefieldnames(options,append);

if nargin == 1
    append = '';
end

cNames = fieldnames(options);
for i = 1:length(cNames)
    eval(['temporaryOptions = options.' cNames{i} ';']);
    if isa(temporaryOptions,'struct')
        cNames = [cNames;recursivefieldnames(temporaryOptions,[cNames{i}])];
    end
end
for i = 1:length(cNames)
    if nargin==1
    else
        cNames{i} = [append '.' cNames{i}];
    end
end




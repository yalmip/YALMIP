function tests = test_scs_review_regressions
tests = functiontests(localfunctions);

function setupOnce(testCase)
% Exercise the interface without requiring or invoking an external solver.
folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
fid = fopen(fullfile(folder.Folder,'scs_direct.m'),'w');
fprintf(fid,'function [x,y,s,info] = scs_direct(data,cones,param)\n');
fprintf(fid,'x = zeros(length(data.c),1); y = zeros(length(data.b),1); s = y;\n');
fprintf(fid,'info.status = param.test_status;\n');
fclose(fid);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder.Folder));

function test_status_mapping(testCase)
m = model;
statuses = {'Solved','Infeasible','Unbounded','Unbounded/Inaccurate','Solved/Inaccurate','Failed'};
expected = [0 1 2 9 9 9];
for i = 1:numel(statuses)
    m.options.scs.test_status = statuses{i};
    out = callscs(m);
    testCase.verifyEqual(out.problem,expected(i));
end

function test_equality_elimination(testCase)
m = model;
m.options.scs.test_status = 'Solved';
m.options.scs.eliminateequalities = 1;
m.c = [1;0]; m.K.f = 1;
m.F_struc = sparse([-2 1 1;3 0 -1]);
out = callscs(m);
testCase.verifyEqual(out.problem,0);
testCase.verifyEqual(sum(out.Primal),2,'AbsTol',1e-12);
testCase.verifyEqual(out.solverinput.cones.z,0);
testCase.verifySize(out.solverinput.data.A,[1 1]);

function m = model
m.F_struc = sparse([1 -1]); m.c = 1;
m.K.f = 0; m.K.l = 1; m.K.q = 0; m.K.e = 0; m.K.p = 0; m.K.s = 0;
m.options.verbose = 0;
m.options.scs = struct('eliminateequalities',0,'gpu',false);
m.options.savedebug = 0; m.options.showprogress = 0;
m.options.savesolverinput = 1; m.options.savesolveroutput = 0;
m.ub = []; m.lb = []; m.evalMap = {}; m.solver.tag = 'scs-direct';

function results = runxunittests(whichtest,xmlFile)

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin
suite = testsuite(whichtest);
runner = TestRunner.withNoPlugins;
if nargin > 1
    p = XMLPlugin.producingJUnitFormat(xmlFile);
else
    [~,name] = fileparts('dev/tests/operators');
    p = XMLPlugin.producingJUnitFormat([name '.xml');
end
runner.addPlugin(p)
results = runner.run(suite);
table(results)
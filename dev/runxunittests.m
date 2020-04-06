function runxunittests(whichtest,xmlFile)

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin
suite = testsuite(whichtest);
runner = TestRunner.withNoPlugins;
if nargin > 1
    p = XMLPlugin.producingJUnitFormat(xmlFile);
else
    p = XMLPlugin.producingJUnitFormat('xunitresults');
end
runner.addPlugin(p)
results = runner.run(suite);
table(results)
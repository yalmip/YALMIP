function runxunittests(whichtest,flag,xmlFile)

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin
suite = testsuite(whichtest);
runner = TestRunner.withNoPlugins;
p = XMLPlugin.producingJUnitFormat(xmlFile);
runner.addPlugin(p)
results = runner.run(suite);
table(results)
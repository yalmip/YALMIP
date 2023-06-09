function results = runxunittests(whichtest,xmlFile)

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin
suite = testsuite(whichtest);
runner = TestRunner.withTextOutput;
if nargin > 1
    p = XMLPlugin.producingJUnitFormat(xmlFile);
else
    [~,name] = fileparts(whichtest);
    xmlFile = [name '.xml'];
    p = XMLPlugin.producingJUnitFormat(xmlFile);
end
runner.addPlugin(p)
results = runner.run(suite);
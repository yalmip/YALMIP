function results = runxunittests(whichtest,xmlFile)

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.XMLPlugin
import matlab.unittest.plugins.CodeCoveragePlugin
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
table(results)
%func_replace_string(xmlFile, xmlFile,'<testsuites>', '')
%func_replace_string(xmlFile, xmlFile,'</testsuites>', '')

function [] = func_replace_string(InputFile, OutputFile, SearchString, ReplaceString)
%%change data [e.g. initial conditions] in model file
% InputFile - string
% OutputFile - string
% SearchString - string
% ReplaceString - string
% read whole model file data into cell array
fid = fopen(InputFile);
data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true);
fclose(fid);
% modify the cell array
% find the position where changes need to be applied and insert new data
for I = 1:length(data{1})
    tf = strcmp(data{1}{I}, SearchString); % search for this string in the array
    if tf == 1
        data{1}{I} = ReplaceString; % replace with this string
    end
end
% write the modified cell array into the text file
fid = fopen(OutputFile, 'w');
for I = 1:length(data{1})
    fprintf(fid, '%s\n', char(data{1}{I}));
end
fclose(fid);
const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const {summarize} = require('./ci-summary.cjs');

function fixture(t, outcomes, files = {}) {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'yalmip-ci-summary-'));
  t.after(() => {
    for (const name of Object.keys(files)) fs.unlinkSync(path.join(directory, name));
    fs.rmdirSync(directory);
  });
  for (const [name, text] of Object.entries(files)) fs.writeFileSync(path.join(directory, name), text);
  const steps = Object.fromEntries(Object.entries(outcomes).map(([id, outcome]) => [id, {outcome}]));
  return summarize(steps, directory, 'integration-basics', 'R2024a', 'https://github.com/example/artifact');
}

test('licensing failure is primary and does not invent test failures', t => {
  const result = fixture(t, {setup: 'success', preflight: 'failure', build: 'skipped', tests: 'skipped'}, {
    'preflight.log': 'License checkout failed.\nLicense Manager Error 4410\nEntitlementNotFound\n',
  });
  assert.match(result.summary, /MATLAB licensing failed/);
  assert.match(result.summary, /EntitlementNotFound/);
  assert.match(result.summary, /Tests did not run/);
  assert.doesNotMatch(result.annotation, /test-reporter|Test execution failed/);
});

test('startup failure still explains the missing tests when MATLAB writes no log', t => {
  const result = fixture(t, {preflight: 'failure', tests: 'skipped'});
  assert.match(result.annotation, /MATLAB startup or the Optimization Toolbox/);
  assert.match(result.summary, /Tests did not run/);
});

test('build errors are shown safely and remain distinct from licensing', t => {
  const result = fixture(t, {preflight: 'success', build: 'failure', tests: 'skipped'}, {
    'sedumi.log': 'mex compilation failed\n<script>not markup</script>\n::error::not an annotation',
  });
  assert.match(result.summary, /Failed at &quot;Build SeDuMi&quot;/);
  assert.match(result.summary, /mex compilation failed/);
  assert.match(result.summary, /&lt;script&gt;/);
  assert.doesNotMatch(result.annotation, /licensing|not an annotation/);
});

test('test failures link to detailed results and downloadable diagnostics', t => {
  const result = fixture(t, {preflight: 'success', build: 'success', tests: 'failure'}, {
    'tests.xml': '<testsuite tests="1" failures="1"/>',
    'tests.log': 'test_operator_example failed',
  });
  assert.match(result.annotation, /Test execution failed/);
  assert.match(result.summary, /MATLAB test summary/);
  assert.match(result.summary, /test_operator_example/);
  assert.match(result.summary, /https:\/\/github.com\/example\/artifact/);
  assert.doesNotMatch(result.summary, /Tests did not run/);
});

test('successful runs need a report and emit no failure annotation', t => {
  const result = fixture(t, {tests: 'success'}, {'tests.xml': '<testsuite tests="1"/>'});
  assert.equal(result.annotation, '');
  assert.equal(result.invalidResults, false);
  assert.match(result.summary, /Test execution succeeded/);
});

test('licensing text inside completed test output does not misclassify a test failure', t => {
  const result = fixture(t, {tests: 'failure'}, {
    'tests.xml': '<testsuite tests="1" failures="1"/>',
    'tests.log': 'Expected diagnostic: License checkout failed.',
  });
  assert.match(result.annotation, /Test execution failed/);
  assert.doesNotMatch(result.annotation, /MATLAB licensing failed/);
});

test('a missing report after purported success is not silently accepted', t => {
  const result = fixture(t, {tests: 'success'});
  assert.equal(result.invalidResults, true);
  assert.match(result.annotation, /No completed test run/);
});

test('cancellation is not misreported as passing or failing tests', t => {
  const result = fixture(t, {tests: 'cancelled'});
  assert.match(result.summary, /cancelled/);
  assert.equal(result.annotation, '');
  assert.equal(result.invalidResults, false);
});

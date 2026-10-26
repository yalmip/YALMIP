// Infrastructure failures must remain visible even when MATLAB cannot write XML.
const fs = require('node:fs');
const path = require('node:path');

const stages = [
  ['checkout', 'Check out YALMIP'],
  ['diagnostics', 'Prepare diagnostics'],
  ['reporting', 'Check CI reporting'],
  ['sedumi', 'Check out SeDuMi'],
  ['revision', 'Record SeDuMi revision'],
  ['setup', 'Install MATLAB and Optimization Toolbox'],
  ['preflight', 'Check MATLAB startup and toolbox license'],
  ['build', 'Build SeDuMi'],
  ['tests', 'Run tests'],
  ['artifacts', 'Upload logs and test results'],
];

function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, c => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  })[c]);
}

function escapeAnnotation(text) {
  return text.replace(/%/g, '%25').replace(/\r/g, '%0D').replace(/\n/g, '%0A');
}

function summarize(steps, directory, group, release, artifactUrl) {
  const outcome = id => steps[id]?.outcome || 'not run';
  const failed = stages.find(([id]) => outcome(id) === 'failure');
  const cancelled = stages.some(([id]) => outcome(id) === 'cancelled');
  const xmlExists = fs.existsSync(path.join(directory, 'tests.xml'));
  let message;
  let excerpt = '';
  let invalidResults = false;

  if (failed) {
    const [id, label] = failed;
    const logfile = {preflight: 'preflight.log', build: 'sedumi.log', tests: 'tests.log'}[id];
    if (logfile && fs.existsSync(path.join(directory, logfile))) {
      const lines = fs.readFileSync(path.join(directory, logfile), 'utf8').split(/\r?\n/);
      const licenseLine = lines.findIndex(line => /license checkout failed|license manager error|licensing error/i.test(line));
      excerpt = (licenseLine >= 0 ? lines.slice(licenseLine, licenseLine + 18) : lines.slice(-40)).join('\n');
    }
    if ((id !== 'tests' || !xmlExists) && /license checkout failed|license manager error|licensing error/i.test(excerpt)) {
      message = `MATLAB licensing failed during "${label}". See the licensing error below and that step's log.`;
    } else if (id === 'preflight') {
      message = 'MATLAB startup or the Optimization Toolbox license check failed. Open "Check MATLAB startup and toolbox license" for the launcher error.';
    } else if (id === 'tests') {
      message = xmlExists
        ? 'Test execution failed. Expand the MATLAB test summary to see failing tests and their diagnostics.'
        : 'Test execution failed before a JUnit report was written. See "Run tests" and tests.log for the startup or test-discovery error.';
    } else {
      message = `Failed at "${label}". Open that step for the original error.`;
    }
    if (outcome('tests') === 'skipped' || outcome('tests') === 'not run') {
      message += ' Tests did not run.';
    }
  } else if (cancelled) {
    message = 'The job was cancelled; results may be incomplete.';
  } else if (outcome('tests') === 'success' && xmlExists) {
    message = 'Test execution succeeded. Detailed results are in the MATLAB test summary.';
  } else {
    message = 'No completed test run with a JUnit report was recorded. Inspect the stage outcomes below.';
    invalidResults = true;
  }

  const lines = [
    `## ${escapeHtml(group)} — ${escapeHtml(release)}`,
    '',
    `**${escapeHtml(message)}**`,
    '',
    '| Stage | Outcome |',
    '| --- | --- |',
    ...stages.map(([id, label]) => `| ${label} | ${escapeHtml(outcome(id))} |`),
    '',
  ];
  if (excerpt) lines.push('### Error context', '', `<pre>${escapeHtml(excerpt.slice(-8000))}</pre>`, '');
  if (artifactUrl) lines.push(`[Download logs, environment details and JUnit results](${artifactUrl}) (retained for 14 days).`, '');
  if (!xmlExists && outcome('tests') !== 'success') {
    lines.push('JUnit results are unavailable because no test report was produced.', '');
  }
  const annotation = failed || invalidResults ? `::error title=CI outcome::${escapeAnnotation(message)}` : '';
  return {summary: lines.join('\n'), annotation, invalidResults};
}

if (require.main === module) {
  const result = summarize(
    JSON.parse(process.env.CI_STEPS || '{}'), 'ci-results',
    process.env.CI_GROUP || 'MATLAB tests', process.env.MATLAB_RELEASE || '',
    process.env.CI_ARTIFACT_URL || '',
  );
  if (process.env.GITHUB_STEP_SUMMARY) fs.appendFileSync(process.env.GITHUB_STEP_SUMMARY, result.summary);
  else console.log(result.summary);
  if (result.annotation) console.log(result.annotation);
  if (result.invalidResults) process.exitCode = 1;
}

module.exports = {summarize};

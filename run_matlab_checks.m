function run_matlab_checks(label)
%RUN_MATLAB_CHECKS The simulator's checks and demo run, in MATLAB.
%   run_matlab_checks(label) runs test_sim, check_saved_scores and a seeded
%   runSimOptimisation, logging to results/matlab_<label>.log and saving the
%   run to results/out_<label>.mat. label defaults to 'run'.
%
%   Start MATLAB in this folder, not inside HPLCMethodOptimisationGUI: the
%   repo's startup.m is the TS-EMO startup script, which adds the TS-EMO
%   folders to the path and warns about a missing Test_functions folder.
%   The path is set by hand rather than with genpath, which would add
%   duplicate copies (Direct/ etc. under Benchmarking/) and put
%   PeakFinder/findpeaks.m ahead of the Signal Toolbox one.

if nargin < 1, label = 'run'; end
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, 'HPLCMethodOptimisationGUI');
addpath(repo, fullfile(repo, 'PeakFinder'), fullfile(repo, 'BayesianOptimiser'), ...
    fullfile(repo, 'BayesianOptimiser', 'GODLIKE'), fullfile(here, 'sim'));
res = fullfile(here, 'results');
if ~exist(res, 'dir'), mkdir(res); end
logFile = fullfile(res, ['matlab_' label '.log']);
if exist(logFile, 'file'), delete(logFile); end
diary(logFile);
cleanup = onCleanup(@() diary('off'));
fprintf('%s  MATLAB %s\n\n', char(datetime('now')), version);

stage('test_sim', @() test_sim(repo));
stage('check_saved_scores', @() check_saved_scores(repo));

rng(1);
out = [];
stage('runSimOptimisation', @() runIt());
if isempty(out), return; end
save(fullfile(res, ['out_' label '.mat']), 'out');
summarise(out);

    function runIt()
        tic;
        out = runSimOptimisation();
        fprintf('runSimOptimisation: %d runs in %.0f s, usedBO = %d\n', numel(out.y), toc, out.usedBO);
    end
end

function stage(name, fn)
fprintf('==== %s\n', name);
try
    fn();
catch ME
    fprintf('%s FAILED: %s\n', name, ME.message);
    for k = 1:min(3, numel(ME.stack))
        fprintf('    at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end
fprintf('\n');
end

function summarise(out)
T = out.truth;
fprintf('==== summary\n');
fprintf('run  GradTime InitialB IsoHold  objective  repoRs  trueRs  nRes  last\n');
for i = 1:numel(out.y)
    fprintf('%3d  %8.2f %8.1f %7.2f  %9.4f  %6.2f  %6.2f  %4d  %5.2f\n', i, out.Z(i, :), ...
        out.y(i), out.repoRsCrit(i), T.RsCrit(i), T.nResolved(i), T.lastPeak(i));
end
b = out.best.byObjective; t = out.best.byTruth;
fprintf('objective''s pick: run %d (true Rs %.2f, %d resolved, last peak %.2f min)\n', ...
    b, T.RsCrit(b), T.nResolved(b), T.lastPeak(b));
fprintf('truth''s pick:     run %d (true Rs %.2f, %d resolved, last peak %.2f min)\n', ...
    t, T.RsCrit(t), T.nResolved(t), T.lastPeak(t));
merged = find(T.RsCrit < 0.5 & out.repoRsCrit > 1);
fprintf('merged but looks resolved (true Rs < 0.5, repo Rs > 1): %d run(s) %s\n', ...
    numel(merged), mat2str(merged'));

% Peaks picked within 0.05 min of t0, split by whether a compound eluted
% there. simHPLC is deterministic, so re-running each method recovers its truth.
sample = simSample();
front = 0; early = 0;
for i = 1:numel(out.y)
    R = simHPLC(out.X(i, :), sample, struct('seed', i));
    P = out.peaksData{i, 5}{1};
    near = P(P(:, 2) > 0 & abs(P(:, 2) - R.t0) < 0.05, 2);
    tR = R.truth.tR(R.truth.eluted);
    for k = 1:numel(near)
        if any(abs(tR - near(k)) < 0.02), early = early + 1; else, front = front + 1; end
    end
end
fprintf('peaks within 0.05 min of t0: %d matching a compound, %d matching none (front)\n', early, front);
end

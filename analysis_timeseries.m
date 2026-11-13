function analysis_timeseries()
%ANALYSIS_TIMESERIES Time-series checks on the simulator's invented output.
%   Three checks, all on simulated data: the noise in one trace, peak
%   measurements repeated over noise seeds, and the objective's history
%   in the Bayesian run saved in results/out_front.mat. Nothing here comes
%   from an instrument.
%
%   Reads results/out_front.mat, so it describes the front-fixed simulator
%   (opts.front = 'measured', commit 28197ad), not the baseline import.
%   Writes nothing except the console log, which is saved as
%   results/timeseries_test.log.
%
%   Provenance and limits, so the numbers are not over-read:
%   A. The noise check tests the pseudo-noise generator in simHPLC
%      (portableNoise): its mean, sd, tail fraction and lag-1/lag-5
%      autocorrelation sit where a Gaussian white process would put them.
%      That is a self-consistency check on the generator. It is not a check
%      of real chromatograph noise, and it says nothing about instruments.
%   B. The replicate spread (sd of retention time about 0.00001 to 0.00004
%      min over 30 seeds) is not evidence of real column repeatability. The
%      noise is deterministic per seed and no run-to-run drift is modelled
%      (see README, Limits). Do not quote it beside a real-instrument
%      repeatability spec without that caveat.
%   C. The objective history is one seeded Bayesian run (seed 1, the same
%      run rerun on the front-fixed simulator). The lag-1 check looks for
%      serial dependence only; it is not a test of the optimiser.

here = fileparts(mfilename('fullpath'));
r = fullfile(here, 'HPLCMethodOptimisationGUI');
addpath(r, fullfile(r, 'PeakFinder'), fullfile(here, 'sim'));
S = simSample();
X = [8 20 95 30 0 0 0 0 2 1.5 1 2 1];
fprintf('Simulated method: GradTime 8 min, %%B 20 to 95, 30 C, 1.5 mL/min, 2 uL\n');
fprintf('Invented sample: %d compounds (simSample)\n\n', numel(S.lnkw));

%% A. Noise in the trace
clean = simHPLC(X, S, struct('noise', 0, 'drift', 0, 'front', 'none'));
z = []; r1 = []; r5 = [];
for s = 1:10
    R = simHPLC(X, S, struct('seed', s, 'drift', 0, 'front', 'none'));
    e = (R.intensity - clean.intensity) / 0.05;   % noise term in units of its sigma
    z = [z; e]; %#ok<AGROW>
    r1(end + 1) = acf(e, 1); %#ok<AGROW>
    r5(end + 1) = acf(e, 5); %#ok<AGROW>
end
N = numel(z) / 10;
fprintf('A. Trace noise (10 seeds, %d points each, residual / 0.05 mAU)\n', N);
fprintf('   mean %.3f, sd %.3f (ideal 0 and 1)\n', mean(z), std(z));
fprintf('   |z| > 3 at %.2f %% of points (Gaussian: 0.27 %%)\n', 100 * mean(abs(z) > 3));
fprintf('   lag-1 autocorrelation %.3f, lag-5 %.3f (95 %% band +/- %.3f)\n\n', ...
    mean(r1), mean(r5), 1.96 / sqrt(N));

%% B. Peak measurements repeated over noise seeds
truth = simHPLC(X, S, struct('noise', 0, 'drift', 0));
tt = truth.truth.tR(truth.truth.eluted);
nC = numel(tt);
nRep = 30;
tR = nan(nRep, nC); H = tR; Ar = tR;
for k = 1:nRep
    R = simHPLC(X, S, struct('seed', k));
    P = simPeaksRow(R); P = P{5}{1};
    P = P(P(:, 2) > 0, :);
    if isempty(P), continue; end
    for c = 1:nC
        [d, j] = min(abs(P(:, 2) - tt(c)));
        if d < 0.05
            tR(k, c) = P(j, 2); H(k, c) = P(j, 3); Ar(k, c) = P(j, 5);
        end
    end
end
fprintf('B. Peaks over %d noise seeds (noise is the only run-to-run term)\n', nRep);
fprintf('   compound  true tR  mean tR  sd tR (min)  found  height CV  area CV\n');
for c = 1:nC
    ok = ~isnan(tR(:, c));
    fprintf('   %8d  %7.4f  %7.4f  %11.5f  %5d  %8.2f %%  %6.2f %%\n', c, tt(c), ...
        mean(tR(ok, c)), std(tR(ok, c)), nnz(ok), ...
        100 * std(H(ok, c)) / mean(H(ok, c)), 100 * std(Ar(ok, c)) / mean(Ar(ok, c)));
end
fprintf('\n');

%% C. Objective history in the Bayesian run
Sf = load(fullfile(here, 'results', 'out_front.mat'));
out = Sf.out;
y = out.y(:);
best = cummin(y);
[~, iBest] = min(y);
fprintf('C. Objective over the %d Bayesian runs (results/out_front.mat)\n', numel(y));
fprintf('   lag-1 autocorrelation of y %.3f (95 %% band +/- %.3f)\n', acf(y, 1), 1.96 / sqrt(numel(y)));
fprintf('   running best improved at %d of %d steps\n', nnz(diff(best) < 0), numel(y) - 1);
fprintf('   best at run %d: objective %.4f, true Rs %.2f\n', iBest, y(iBest), out.truth.RsCrit(iBest));
end

function r = acf(x, lag)
x = x(:) - mean(x(:));
r = sum(x(1:end - lag) .* x(1 + lag:end)) / sum(x .^ 2);
end

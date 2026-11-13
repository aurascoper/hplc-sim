function export_timeseries_csv()
%EXPORT_TIMESERIES_CSV Long-format CSVs of the simulator's invented series, for the dashboard.
%   Reads results/out_front.mat (the front-fixed simulator, commit 28197ad)
%   and writes results/ts_*.csv. These CSVs are a direct export of that run,
%   not an independent measurement: every trace, noise series, replicate and
%   history value is recomputed from the simulator with the same seeds that
%   runSimOptimisation uses, so they match out_front.mat exactly. Nothing is
%   timestamped and no path is written into the files, so re-exporting
%   gives byte-identical CSVs unless the simulator or the run changes.
%   The dashboard page reads these files; treat its numbers as the
%   simulator's output, not as separate data.

here = fileparts(mfilename('fullpath'));
r = fullfile(here, 'HPLCMethodOptimisationGUI');
addpath(r, fullfile(r, 'PeakFinder'), fullfile(here, 'sim'));
res = fullfile(here, 'results');
S = simSample();
Sf = load(fullfile(res, 'out_front.mat'));
out = Sf.out;
names = {'A','B','C','D','E','F','G','H'};

% traces and %B programs for four runs of the Bayesian history
picks = [1 12 16 27];
fid = fopen(fullfile(res, 'ts_traces.csv'), 'w'); fprintf(fid, 'run,time_min,intensity_mAU\n');
fp = fopen(fullfile(res, 'ts_program.csv'), 'w'); fprintf(fp, 'run,time_min,pctB\n');
for i = picks
    R = simHPLC(out.X(i, :), S, struct('seed', i, 'front', 'measured'));
    idx = 1:2:numel(R.time);
    for j = idx
        fprintf(fid, '%d,%.4f,%.4f\n', i, R.time(j), R.intensity(j));
    end
    for j = 1:size(R.prog, 1)
        fprintf(fp, '%d,%.4f,%.2f\n', i, R.prog(j, 1), R.prog(j, 2));
    end
end
fclose(fid); fclose(fp);

% noise: residual of the trace from its noise-free version, in units of its sigma
X1 = [8 20 95 30 0 0 0 0 2 1.5 1 2 1];
clean = simHPLC(X1, S, struct('noise', 0, 'drift', 0, 'front', 'none'));
z = []; acfs = nan(10, 41);
for s = 1:10
    R = simHPLC(X1, S, struct('seed', s, 'drift', 0, 'front', 'none'));
    e = (R.intensity - clean.intensity) / 0.05;
    z = [z; e]; %#ok<AGROW>
    if s == 1, z1 = e; end
    for lag = 0:40
        acfs(s, lag + 1) = acf(e, lag);
    end
end
N = numel(z);
edges = -4:0.25:4;
counts = histcounts(z, edges);
fid = fopen(fullfile(res, 'ts_noise_hist.csv'), 'w'); fprintf(fid, 'bin_low,bin_high,observed,expected\n');
for b = 1:numel(counts)
    p = 0.5 * (erf(edges(b+1) / sqrt(2)) - erf(edges(b) / sqrt(2)));
    fprintf(fid, '%.2f,%.2f,%d,%.2f\n', edges(b), edges(b+1), counts(b), N * p);
end
fclose(fid);
fid = fopen(fullfile(res, 'ts_noise_acf.csv'), 'w'); fprintf(fid, 'lag,acf,band\n');
for lag = 0:40
    fprintf(fid, '%d,%.4f,%.4f\n', lag, mean(acfs(:, lag + 1)), 1.96 / sqrt(numel(z1)));
end
fclose(fid);
fid = fopen(fullfile(res, 'ts_noise_seed1.csv'), 'w'); fprintf(fid, 'sample,z\n');
for j = 1:numel(z1)
    fprintf(fid, '%d,%.4f\n', j, z1(j));
end
fclose(fid);

% replicates over 30 noise seeds (same method as analysis_timeseries.m)
truth = simHPLC(X1, S, struct('noise', 0, 'drift', 0));
tt = truth.truth.tR(truth.truth.eluted);
nC = numel(tt); nRep = 30;
fid = fopen(fullfile(res, 'ts_replicates.csv'), 'w'); fprintf(fid, 'seed,compound,tR_min,height_mAU,area,true_tR_min\n');
for k = 1:nRep
    R = simHPLC(X1, S, struct('seed', k));
    P = simPeaksRow(R); P = P{5}{1};
    P = P(P(:, 2) > 0, :);
    for c = 1:nC
        if isempty(P), break; end
        [d, j] = min(abs(P(:, 2) - tt(c)));
        if d < 0.05
            fprintf(fid, '%d,%s,%.5f,%.3f,%.4f,%.5f\n', k, names{c}, P(j, 2), P(j, 3), P(j, 5), tt(c));
        end
    end
end
fclose(fid);

% objective history of the Bayesian run
y = out.y(:);
best = cummin(y);
fid = fopen(fullfile(res, 'ts_history.csv'), 'w');
fprintf(fid, 'run,GradTime_min,InitialB_pct,IsoHold_min,objective,best_so_far,true_RsCrit,repo_RsCrit,n_resolved,last_peak_min\n');
for i = 1:numel(y)
    fprintf(fid, '%d,%.3f,%.2f,%.3f,%.5f,%.5f,%.4f,%.4f,%d,%.3f\n', i, out.Z(i,1), out.Z(i,2), ...
        out.Z(i,3), y(i), best(i), out.truth.RsCrit(i), out.repoRsCrit(i), out.truth.nResolved(i), out.truth.lastPeak(i));
end
fclose(fid);
fprintf('wrote results/ts_*.csv from results/out_front.mat\n');
end

function c = acf(x, lag)
x = x(:) - mean(x(:));
c = sum(x(1:end - lag) .* x(1 + lag:end)) / sum(x .^ 2);
end

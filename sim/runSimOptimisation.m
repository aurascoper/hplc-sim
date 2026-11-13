function out = runSimOptimisation(varargin)
%RUNSIMOPTIMISATION Closed-loop method optimisation against the simulator.
%   out = runSimOptimisation() runs the repo's loop with simHPLC in place of
%   ChemStation: an initial design, then one suggested method per step,
%   each scored by the repo's own objective function on the repo's own
%   peak picking. It records the simulator's ground truth beside every
%   score, which a real instrument cannot give you.
%
%   Name-value options (defaults):
%     'objective'  @WeightRespcmtmdHPLC018 (what Optimisations 1-3 used)
%     'nInit'      7      initial experiments (Optimisation 1 used 7)
%     'nOpt'       25     optimisation steps after that
%     'sample'     simSample()
%     'seed'       1      seeds the initial design, the fallback search and
%                         the global random stream that BayesianOptimiser
%                         and GODLIKE draw from, so a run repeats exactly;
%                         the caller's random state is restored on return
%     'useBO'      true   use BayesianOptimiser (needs fitrgp, from the
%                         Statistics and Machine Learning Toolbox); with
%                         false, or without fitrgp, random search instead
%     'front'      'measured'  simHPLC's opts.front; 'none' removes the
%                         solvent front, to see what voidTime 0 lets through
%   Design space, as in Optimisation 1: GradTime 1-10 min, initial %B 5-60,
%   iso hold 0-10 min; final %B 95, 30 C, 2 uL, 1.5 mL/min fixed.
%
%   out.X, out.y          methods run and their objective values
%   out.truth             per run: RsCrit (true, USP), nResolved, lastPeak
%   out.repoRsCrit        the repo formula's critical resolution per run
%   out.best              index of the run the objective ranks best, and the
%                         run the truth ranks best (fewest unresolved peaks,
%                         then shortest)
%
%   The objective is called the way the app calls it, with voidTime 0 and
%   skew 1.1. With voidTime 0 the void-time filter removes nothing, so the
%   measured solvent front counts as a peak, as it did in all 420 saved runs.

p = struct('objective', @WeightRespcmtmdHPLC018, 'nInit', 7, 'nOpt', 25, ...
    'sample', [], 'seed', 1, 'useBO', true, 'front', 'measured');
for i = 1:2:numel(varargin)
    p.(varargin{i}) = varargin{i + 1};
end
if isempty(p.sample), p.sample = simSample(); end
callerRng = rng;
restoreRng = onCleanup(@() rng(callerRng)); %#ok<NASGU>
rng(p.seed);
lb = [1 5 0];
ub = [10 60 10];
toMethod = @(z) [z(1) z(2) 95 30 z(3) 0 0 0 2 1.5 1 2 1];

useBO = p.useBO && exist('fitrgp', 'file') == 2 && exist('BayesianOptimiser', 'file') == 2;
if p.useBO && ~useBO
    fprintf('runSimOptimisation: fitrgp or BayesianOptimiser not found; using random search\n');
end

Z = designLHS(p.nInit, lb, ub, p.seed);
peaksData = cell(0, 5);
truth = struct('RsCrit', {}, 'nResolved', {}, 'lastPeak', {}, 'nEluted', {});
for i = 1:size(Z, 1)
    [peaksData, truth] = runOne(toMethod(Z(i, :)), p, peaksData, truth, i);
end
y = scoreAll(p.objective, peaksData);

for s = 1:p.nOpt
    if useBO
        opt = BayesianOptimiser('ardmatern52', 'AEI', [lb; ub], Z, y);
        z = opt.suggest();
    else
        z = designLHS(1, lb, ub, p.seed + 1000 + s);
    end
    Z(end + 1, :) = z; %#ok<AGROW>
    [peaksData, truth] = runOne(toMethod(z), p, peaksData, truth, size(Z, 1));
    y = scoreAll(p.objective, peaksData);
end

repoRs = zeros(size(Z, 1), 1);
for i = 1:size(Z, 1)
    P = peaksData{i, 5}{1};
    P = P(P(:, 2) > 0, :);
    repoRs(i) = repoCrit(P);
end
T = struct('RsCrit', [truth.RsCrit]', 'nResolved', [truth.nResolved]', ...
    'lastPeak', [truth.lastPeak]', 'nEluted', [truth.nEluted]');
[~, bestObj] = min(y);
key = [-T.nResolved, T.lastPeak];
[~, ord] = sortrows(key);
out = struct('Z', Z, 'X', cell2mat(arrayfun(@(i) toMethod(Z(i, :)), (1:size(Z, 1))', 'UniformOutput', false)), ...
    'y', y, 'truth', T, 'repoRsCrit', repoRs, 'usedBO', useBO, ...
    'best', struct('byObjective', bestObj, 'byTruth', ord(1)), 'peaksData', {peaksData});
end

function [peaksData, truth] = runOne(X, p, peaksData, truth, i)
R = simHPLC(X, p.sample, struct('seed', i, 'front', p.front));
peaksData(end + 1, :) = simPeaksRow(R);
truth(end + 1) = struct('RsCrit', R.truth.RsCrit, 'nResolved', R.truth.nResolved, ...
    'lastPeak', R.truth.lastPeak, 'nEluted', R.truth.nEluted);
end

function y = scoreAll(fn, peaksData)
% The app recomputes every score when a run is added (calcAllRespFunc),
% because the combined objectives normalise by the longest run so far.
n = size(peaksData, 1);
y = zeros(n, 1);
for i = 1:n
    y(i) = feval(fn, peaksData, 0, 1.1, i, 1);
end
end

function rs = repoCrit(P)
if size(P, 1) < 2, rs = NaN; return; end
rt = P(:, 2); w = P(:, 4) * 2.355;
rs = inf;
for a = 1:numel(rt)
    for b = a + 1:numel(rt)
        rs = min(rs, 2 * abs(rt(b) - rt(a)) / (w(a) + w(b)));
    end
end
end

function Z = designLHS(n, lb, ub, seed)
% A Latin hypercube that does not depend on lhsdesign or the host RNG.
d = numel(lb);
Z = zeros(n, d);
for j = 1:d
    u = mod(sin((1:n)' * (12.9898 + j) + seed * 78.233) * 43758.5453, 1);
    [~, perm] = sort(mod(sin((1:n)' * (3.1 + 7 * j) + seed * 5.77) * 9631.27, 1));
    Z(:, j) = lb(j) + (ub(j) - lb(j)) * ((perm - 1) + u) / n;
end
end

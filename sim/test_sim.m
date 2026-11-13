function test_sim(repoRoot)
%TEST_SIM Checks for the simulated chromatograph. MATLAB or Octave.
%   test_sim(repoRoot) needs the HPLCMethodOptimisationGUI checkout (for
%   PeakFinder/autofindpeaks and the objective functions) at repoRoot.
%
%   1  isocratic retention matches t0*(1+k) and t0 = Vm/F
%   2  linear-gradient retention matches the closed-form LSS result
%   3  methodTimetable reproduces RunCond's three timetable branches
%   4  autofindpeaks recovers the simulated positions and FWHM
%   5  the repo's resolution formula reads 0.721 of the true resolution,
%      because it multiplies FWHM by 2.355 where base width is 1.699 FWHM
%   6  temperature lowers retention when dH < 0, and Column/OMPhase/pH
%      change nothing (they are declared ignored)

here = fileparts(mfilename('fullpath'));
addpath(here);
if nargin >= 1 && ~isempty(repoRoot)
    addpath(repoRoot);
    addpath(fullfile(repoRoot, 'PeakFinder'));
end
quiet = struct('noise', 0, 'drift', 0);
n = 0;

% 1 isocratic: B fixed at 50 % for the whole run (iso hold longer than the run)
pair = simSample('pair');
X = method(10, 50, 50, 30, 30, 0, 0, 0, 2, 1.0);
R = simHPLC(X, pair, quiet);
k = exp(pair.lnkw - pair.S * 0.5);
expected = R.t0 * (1 + k);
check(abs(R.t0 - 0.5) < 1e-12, 't0 = Vm/F');
check(all(abs(R.truth.tR - expected) ./ expected < 0.002), sprintf( ...
    'isocratic tR %s vs t0(1+k) %s', mat2str(R.truth.tR', 5), mat2str(expected', 5)));
n = n + 2;

% 2 linear gradient from 5 % to 95 % over 20 min, no dwell volume
o = quiet; o.Vd = 0;
X = method(20, 5, 95, 30, 0, 0, 0, 0, 2, 1.0);
R = simHPLC(X, pair, o);
beta = 0.90 / 20;
k0 = exp(pair.lnkw - pair.S * 0.05);
b = pair.S * beta * R.t0;
lss = (R.t0 ./ b) .* log(1 + b .* k0) + R.t0;
check(all(abs(R.truth.tR - lss) ./ lss < 0.01), sprintf( ...
    'gradient tR %s vs LSS %s', mat2str(R.truth.tR', 5), mat2str(lss', 5)));
n = n + 1;

% 3 timetables, against RunCond.m's rows written out by hand
[p, m] = methodTimetable(method(4, 20, 80, 30, 0, 0, 0, 0, 2, 1));
check(isequal(p, [0 20; 4 80; 4.1 95; 10 95]) && m == 10, 'branch 1 timetable');
[p, m] = methodTimetable(method(4, 20, 80, 30, 2, 0, 0, 0, 2, 1));
check(isequal(p, [0 20; 2 20; 6 80; 6.1 95; 11 95]) && m == 12, 'branch 2 timetable (stop +6, hold to +5, as RunCond)');
[p, m] = methodTimetable(method(4, 20, 80, 30, 2, 1, 3, 50, 2, 1));
check(isequal(p, [0 20; 2 20; 6 50; 7 50; 10 80; 10.1 95; 16 95]) && m == 16, 'branch 3 timetable');
n = n + 3;

% 4 peak picking on a clean trace
if exist('autofindpeaks', 'file') == 2
    X = method(8, 20, 95, 30, 0, 0, 0, 0, 2, 1.5);
    R = simHPLC(X, simSample(), struct('noise', 0.02));
    row = simPeaksRow(R);
    P = row{5}{1};
    P = P(P(:, 2) > R.t0 + 0.1, :);
    tr = sort(R.truth.tR(R.truth.eluted));
    resolved = R.truth.RsAdjacent;
    for i = 1:numel(tr)
        isolated = (i == 1 || resolved(i - 1) > 1.5) && (i == numel(tr) || resolved(i) > 1.5);
        if ~isolated, continue; end
        [d, j] = min(abs(P(:, 2) - tr(i)));
        check(d < 0.01, sprintf('peak at %.3f found within 0.01 min (nearest %.3f)', tr(i), P(j, 2)));
        f = R.truth.fwhm(R.truth.tR == tr(i));
        check(abs(P(j, 4) - f) / f < 0.15, sprintf('FWHM %.4f vs true %.4f', P(j, 4), f));
        n = n + 2;
    end
else
    fprintf('test_sim: autofindpeaks not on the path; peak-picking checks skipped\n');
end

% 5 the resolution factor
sig = [0.02 0.02];
t = [2.0 2.2];
fwhm = 2 * sqrt(2 * log(2)) * sig;
trueRs = 2 * (t(2) - t(1)) / (4 * sig(1) + 4 * sig(2));
repoRs = 2 * (t(2) - t(1)) / (2.355 * fwhm(1) + 2.355 * fwhm(2));
check(abs(repoRs / trueRs - 0.7213) < 1e-3, sprintf('repo Rs / true Rs = %.4f', repoRs / trueRs));
n = n + 1;

% 6 temperature and the ignored variables
Xa = method(8, 20, 95, 30, 0, 0, 0, 0, 2, 1.5);
Xb = Xa; Xb(4) = 50;
Ra = simHPLC(Xa, simSample(), quiet);
Rb = simHPLC(Xb, simSample(), quiet);
check(all(Rb.truth.tR(Rb.truth.eluted) < Ra.truth.tR(Rb.truth.eluted)), 'hotter elutes earlier');
Xc = Xa; Xc(11:13) = [3 1 7];
Rc = simHPLC(Xc, simSample(), quiet);
check(isequal(Rc.intensity, Ra.intensity), 'Column, OMPhase and pH change nothing');
n = n + 2;

fprintf('test_sim: %d checks passed\n', n);
end

function X = method(G, B0, B1, T, iso, iso2, G2, sw, inj, F)
X = [G B0 B1 T iso iso2 G2 sw inj F 1 2 1];
end

function check(ok, label)
if ~ok
    error('test_sim:fail', '%s', label);
end
end

function R = simHPLC(X, sample, opts)
%SIMHPLC Simulated chromatograph standing in for ChemStation.
%   R = simHPLC(X, sample, opts) runs one method (a row of RunCond's 13
%   variables) on a sample from simSample and returns
%     R.time, R.intensity   the detector trace (min, mAU), one channel
%     R.prog, R.maxTime     the %B program from methodTimetable
%     R.t0, R.tD            column void time and dwell time (min)
%     R.truth               per compound: tR, sigma, fwhm, area, eluted,
%                           kElution; and the true adjacent resolutions
%
%   The physics, all of it standard and all of it simplified:
%     retention    ln k = lnkw(T) - S*phi, with van 't Hoff for temperature:
%                  lnkw(T) = lnkw - dH/R*(1/T - 1/T0), T in kelvin
%     migration    dx/dt = 1/(t0*(1+k)) for column fraction x, where k is
%                  taken at the composition that reached that point of the
%                  column: phi(t - tD - x*t0). Integrated numerically, so
%                  multi-step gradients need no special case.
%     band width   sigma = t0*(1+k_elution)/sqrt(N). No gradient band
%                  compression, so gradient peaks come out wider than real
%                  ones: the simulator is conservative about resolution.
%     plates       N(F) = N0*H(F0)/H(F), H = A + B/F + C*F, a van Deemter
%                  shape with invented A, B, C.
%     detector     exponentially modified Gaussian-like peaks with slight
%                  tailing, a solvent-front blip at t0, baseline drift with
%                  %B, white noise.
%   A compound still on the column at the stop time is reported with
%   eluted = false and contributes nothing to the trace.
%
%   Column, OMPhase and pH (variables 11-13) are ignored: refusing them
%   would block the app's default design, and modelling them would invent
%   more than this simulator already does. Change them and nothing moves.
%
%   opts fields (defaults): Vm 0.5 mL, Vd 0.2 mL, N0 4000 at F0 1 mL/min,
%   vanDeemter [1 0.4 0.6], rate 150 points/min, noise 0.05 mAU,
%   drift 0.02 mAU per %B, seed 1, dt 0.002 min.

if nargin < 2 || isempty(sample), sample = simSample(); end
if nargin < 3, opts = struct(); end
o = defaults(opts);
if iscell(X), X = cell2mat(X); end
[prog, maxTime] = methodTimetable(X);
T = X(4); InjVol = X(9); F = X(10);
if ~(F > 0), error('hplcsim:simHPLC:flow', 'Flow rate must be positive.'); end

t0 = o.Vm / F;
tD = o.Vd / F;
Rgas = 8.314462618;
lnkwT = sample.lnkw - sample.dH / Rgas * (1 / (T + 273.15) - 1 / (sample.T0 + 273.15));
H = @(f) o.vanDeemter(1) + o.vanDeemter(2) ./ f + o.vanDeemter(3) .* f;
N = o.N0 * H(o.F0) / H(F);

phiAt = @(t) interp1([-1e9; prog(:, 1); 1e9], [prog(1, 2); prog(:, 2); prog(end, 2)], t, 'linear') / 100;

n = numel(sample.lnkw);
x = zeros(n, 1);
tR = nan(n, 1);
kE = nan(n, 1);
t = 0;
while t < maxTime && any(isnan(tR))
    live = isnan(tR);
    phi = phiAt(t - tD - x(live) * t0);
    k = exp(lnkwT(live) - sample.S(live) .* phi);
    xNew = x;
    xNew(live) = x(live) + o.dt ./ (t0 * (1 + k));
    done = find(live & xNew >= 1);
    for i = done'
        frac = (1 - x(i)) / (xNew(i) - x(i));
        tR(i) = t + frac * o.dt;
        kE(i) = exp(lnkwT(i) - sample.S(i) * phiAt(tR(i) - tD - t0));
    end
    x = xNew;
    t = t + o.dt;
end
eluted = ~isnan(tR) & tR <= maxTime;
sigma = t0 * (1 + kE) / sqrt(N);
area = sample.rf * InjVol;

time = (1 / o.rate:1 / o.rate:maxTime)';
y = zeros(size(time));
tau = 0.1 * sigma;
for i = find(eluted)'
    d = time - tR(i);
    den = 2 * sigma(i)^2 + tau(i) * d;
    shape = zeros(size(d));
    ok = den > 0;
    shape(ok) = exp(-d(ok).^2 ./ den(ok));
    s = trapz(time, shape);
    if s > 0
        y = y + area(i) * shape / s;
    end
end
% Solvent front: a small positive then negative blip at t0.
y = y + 3 * exp(-(time - t0).^2 / (2 * 0.01^2)) - 2 * exp(-(time - t0 - 0.03).^2 / (2 * 0.015^2));
y = y + o.drift * 100 * phiAt(time - tD - t0);
y = y + o.noise * portableNoise(numel(time), o.seed);

order = find(eluted);
[~, srt] = sort(tR(order));
order = order(srt);
RsTrue = nan(max(numel(order) - 1, 0), 1);
for j = 1:numel(order) - 1
    a = order(j); b = order(j + 1);
    RsTrue(j) = 2 * (tR(b) - tR(a)) / (4 * sigma(a) + 4 * sigma(b));
end

R = struct();
R.name = sprintf('SIM %s T%.0f F%.2f G%.1f B%.0f-%.0f', sample.name, T, F, X(1), X(2), X(3));
R.time = time;
R.intensity = y;
R.prog = prog;
R.maxTime = maxTime;
R.t0 = t0;
R.tD = tD;
R.N = N;
R.truth = struct('tR', tR, 'sigma', sigma, 'fwhm', 2 * sqrt(2 * log(2)) * sigma, ...
    'area', area, 'eluted', eluted, 'kElution', kE, 'elutionOrder', order, ...
    'RsAdjacent', RsTrue, 'RsCrit', minOrNaN(RsTrue), 'lastPeak', maxOrNaN(tR(eluted)), ...
    'nEluted', nnz(eluted), 'nResolved', countResolved(RsTrue, nnz(eluted)));
end

function o = defaults(o)
d = struct('Vm', 0.5, 'Vd', 0.2, 'N0', 4000, 'F0', 1, 'vanDeemter', [1 0.4 0.6], ...
    'rate', 150, 'noise', 0.05, 'drift', 0.02, 'seed', 1, 'dt', 0.002);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i}), o.(f{i}) = d.(f{i}); end
end
end

function v = minOrNaN(x)
if isempty(x), v = NaN; else, v = min(x); end
end

function v = maxOrNaN(x)
if isempty(x), v = NaN; else, v = max(x); end
end

function n = countResolved(Rs, nEluted)
% Peaks that stand apart from both neighbours at Rs >= 1.5.
if nEluted == 0, n = 0; return; end
ok = true(nEluted, 1);
bad = find(Rs < 1.5);
ok(bad) = false;
ok(bad + 1) = false;
n = nnz(ok);
end

function z = portableNoise(n, seed)
% Deterministic standard-normal-like noise that does not depend on the
% host's random generator, so MATLAB and Octave draw the same kind of
% sequence. Not for statistics; for a realistic baseline only.
i = (1:n)';
u1 = mod(sin(i * 12.9898 + seed * 78.233) * 43758.5453, 1);
u2 = mod(sin(i * 39.3468 + seed * 11.135) * 24634.6345, 1);
u1 = min(max(u1, 1e-12), 1);
z = sqrt(-2 * log(u1)) .* cos(2 * pi * u2);
end

function S = simSample(name)
%SIMSAMPLE An invented test mixture for the simulated chromatograph.
%   S = simSample() or simSample('mix8') returns eight compounds, each with
%     lnkw  ln k extrapolated to 0 % B at 30 C
%     S     slope of ln k against the B fraction (natural-log basis)
%     dH    retention enthalpy, J/mol (negative: retention falls with T)
%     rf    peak area per microlitre injected, mAU min / uL
%
%   Every value is invented. The set is shaped to be a useful exercise, not
%   to stand for any real sample. Two pairs set the difficulty:
%     4 and 5   swap order near 19 % B (ln k equal where 0.3 = 1.6*phi)
%     6 and 7   merge near 60 % B (0.3 = 0.5*phi), and differ in dH, so temperature moves
%               them relative to each other
%   Retention spans weak to strong, so a short gradient loses peaks and a
%   long one wastes time. Tuned by a coarse scan of the Optimisation 1
%   design space (GradTime 1-10, initial %B 5-60, iso hold 0-10; 30 C,
%   1.5 mL/min): about 8 % of that grid resolves all eight at Rs >= 1.5.
%
%   S = simSample('pair') returns only the critical pair, for tests.

if nargin < 1, name = 'mix8'; end
%        lnkw   S      dH(J/mol)  rf
T = [    2.6    9.0    -9000     0.9
         3.4   10.5   -11000     1.2
         4.3   12.0   -12000     0.8
         5.0   12.0   -14000     1.0
         5.3   13.6   -14000     1.1
         6.1   14.0   -20000     0.7
         6.4   14.5   -12000     0.9
         7.4   15.5   -16000     1.3 ];
switch name
    case 'mix8'
    case 'pair'
        T = T(4:5, :);
    otherwise
        error('hplcsim:simSample:name', 'Unknown sample "%s".', name);
end
S = struct('name', name, 'lnkw', T(:, 1), 'S', T(:, 2), 'dH', T(:, 3), 'rf', T(:, 4), ...
    'T0', 30, 'note', 'invented test mixture; not a real sample');
end

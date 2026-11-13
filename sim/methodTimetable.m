function [prog, maxTime] = methodTimetable(X)
%METHODTIMETABLE The %B program RunCond.m would write for one experiment.
%   [prog, maxTime] = methodTimetable(X) takes one row of the 13 method
%   variables in RunCond's order:
%     1 GradTime  2 BInitial  3 BFinal  4 Temp  5 IsoTime  6 IsoTime2
%     7 GradTime2 8 Switch    9 InjVol 10 FlowRate 11 Column 12 OMPhase 13 pH
%   and returns prog = [time_min, percentB] rows, starting at t = 0, and
%   the run's stop time (cell E1 in RunCond).
%
%   It copies RunCond's three branches rather than tidying them, so that
%   the simulator runs the method the real system would have run:
%     IsoTime == 0 and IsoTime2 == 0   gradient, ramp to 95 %, hold to +6 min
%     IsoTime2 == 0                    hold, gradient, ramp, hold to +5 min
%     otherwise                        hold, gradient to the switch point,
%                                      hold, second gradient, ramp, hold to +6
%   Two RunCond details carried over as written: in the first two branches
%   the stop time is GradTime + IsoTime + 6 even where the 95 % hold ends at
%   +5, and the 0.1 min ramp to 95 % is part of the program.
%   Column, OMPhase and pH (11-13) are not part of the timetable.

if iscell(X)
    X = cell2mat(X);
end
if numel(X) < 10
    error('hplcsim:methodTimetable:input', 'Need at least the first 10 method variables.');
end
G = X(1); B0 = X(2); B1 = X(3); iso = X(5); iso2 = X(6); G2 = X(7); sw = X(8);

if iso == 0 && iso2 == 0
    prog = [0 B0; G B1; G + 0.1 95; G + 6 95];
    maxTime = G + iso + 6;
elseif iso2 == 0
    prog = [0 B0; iso B0; G + iso B1; G + iso + 0.1 95; G + iso + 5 95];
    maxTime = G + iso + 6;
else
    Bsw = (sw / 100) * (B1 - B0) + B0;
    t1 = iso + G; t2 = t1 + iso2; t3 = t2 + G2;
    prog = [0 B0; iso B0; t1 Bsw; t2 Bsw; t3 B1; t3 + 0.1 95; t3 + 6 95];
    maxTime = t3 + 6;
end
if any(diff(prog(:, 1)) < 0)
    error('hplcsim:methodTimetable:order', 'The timetable runs backwards; check the method variables.');
end
end

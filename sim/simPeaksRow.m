function row = simPeaksRow(R, settings)
%SIMPEAKSROW One experiment in the layout the app's objectives read.
%   row = simPeaksRow(R, settings) peak-picks the simulated trace with the
%   repo's own PeakFinder/autofindpeaks and returns a 1x5 cell laid out as
%   the saved app data are (PaperResults/*.mat, data row "6"):
%     {name, time, wavelengths, intensity, {peakTable}}
%   so peaksData{currExp,5}{:, 1} is the peak table
%   [number, position, height, FWHM, area] the objective functions expect.
%
%   settings defaults are the ones saved with Optimisation 1:
%   ampThresh 10, slopeThresh 0.05, smoothWidth 1, peakGroup 5. The
%   single simulated channel is labelled 210 nm, the saved WavelengthSingle.

if nargin < 2, settings = struct(); end
d = struct('ampThresh', 10, 'slopeThresh', 0.05, 'smoothWidth', 1, 'peakGroup', 5, 'wavelength', 210);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(settings, f{i}), settings.(f{i}) = d.(f{i}); end
end
if exist('autofindpeaks', 'file') ~= 2
    error('hplcsim:simPeaksRow:path', 'autofindpeaks is not on the path; addpath the repo''s PeakFinder folder.');
end
P = autofindpeaks(R.time, R.intensity, settings.slopeThresh, settings.ampThresh, ...
    settings.smoothWidth, settings.peakGroup, 3);
% autofindpeaks returns [0 0 0 0 0] when it finds nothing; the saved data
% keep that row, and the objectives treat it as one peak at t = 0, which
% the void-time filter then removes.
row = {R.name, R.time, settings.wavelength, R.intensity, {P}};
end

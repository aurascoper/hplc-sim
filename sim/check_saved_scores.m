function report = check_saved_scores(repoRoot)
%CHECK_SAVED_SCORES Recompute the paper's saved objective scores.
%   report = check_saved_scores(repoRoot) loads each
%   PaperResults/*/Optimisation*.mat, re-runs the repo's five objective
%   functions on the saved peak tables (at the saved single wavelength,
%   with voidTime 0 and skew 1.1, as HPLCapp calls them) and compares
%   with the saved InitalResults table.
%
%   Optimisations 1-3 were saved before the three single objectives gained
%   their log transforms, so for those the old form is also tried:
%     RsCrit2              -(Rs - 2)       now  -log(Rs)
%     Maximisepeaknumber   -n              now  -log(n)
%     Minimisemethodtime   t_last          now  log(t_last)
%   and the report says which form matched.
%
%   Run in Octave on 2026-11-12 from a Python export of the same files
%   (Octave cannot read MATLAB tables): all 420 runs matched to 4e-15,
%   except CombinedFunction1 in Optimisations 4 and 7 (0.013, 0.002). That
%   objective normalises by the longest run so far, so a score saved before
%   the last runs existed differs. This MATLAB version reads the .mat files
%   directly and has not itself been run.

names = {'RsCrit2', 'Maximisepeaknumber', 'Minimisemethodtime', 'CombinedFunction1', 'WeightRespcmtmdHPLC018'};
fns = {@responceRscrit, @responceNumbPeaks, @responceMaxTime, @responceCombFunc1, @WeightRespcmtmdHPLC018};
old = {@(v) -(exp(-v) - 2), @(v) -exp(-v), @(v) exp(v), @(v) v, @(v) v};

files = dir(fullfile(repoRoot, 'PaperResults', '*', 'Optimisation*.mat'));
report = struct('file', {}, 'nRuns', {}, 'maxDevCurrent', {}, 'maxDevOldForm', {});
for f = 1:numel(files)
    S = load(fullfile(files(f).folder, files(f).name));
    d = S.data;
    key = cellfun(@char, d(:, 1), 'UniformOutput', false);
    D = d{strcmp(key, 'Data'), 2};
    res = d{strcmp(key, 'InitalResults'), 2};
    st = d{strcmp(key, 'Settings'), 2};
    n = size(D, 1);
    devNow = zeros(1, 5);
    devOld = zeros(1, 5);
    for e = 1:n
        ind = find(D{e, 3} == st.WavelengthSingle);
        for j = 1:5
            saved = res.(names{j})(e);
            if iscell(saved), saved = saved{1}; end
            v = feval(fns{j}, D, 0, 1.1, e, ind);
            devNow(j) = max(devNow(j), abs(v - saved));
            devOld(j) = max(devOld(j), abs(old{j}(v) - saved));
        end
    end
    report(end + 1) = struct('file', files(f).name, 'nRuns', n, ...
        'maxDevCurrent', devNow, 'maxDevOldForm', devOld); %#ok<AGROW>
    fprintf('%-20s %2d runs  current %s  old form %s\n', files(f).name, n, ...
        mat2str(devNow, 3), mat2str(devOld, 3));
end
end

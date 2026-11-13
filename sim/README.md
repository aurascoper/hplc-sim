# Simulated chromatograph for HPLCMethodOptimisationGUI

A stand-in for ChemStation, so the Bourne Group's closed-loop optimiser can run at home with no instrument. The method writer, peak picking and objective functions stay the repo's own. Only the instrument is replaced.

Put this folder beside the repo checkout, then in MATLAB:

```matlab
addpath(genpath('HPLCMethodOptimisationGUI'))   % repo: PeakFinder, objectives, BayesianOptimiser
addpath('sim')
test_sim('HPLCMethodOptimisationGUI')            % 21 checks
check_saved_scores('HPLCMethodOptimisationGUI')  % the paper's 420 saved runs
out = runSimOptimisation();                      % 7 initial runs + 25 Bayesian steps
```

`runSimOptimisation` needs `fitrgp` (Statistics and Machine Learning Toolbox) for the repo's `BayesianOptimiser`. Without it, it falls back to random search and says so. In Octave, also `addpath('sim/octave')` for a `fullfact` stand-in.

| File | Does |
|---|---|
| `methodTimetable.m` | The %B program `RunCond.m` would write, copied branch for branch, including the stop time running past the 95 % hold in the second branch. |
| `simSample.m` | An invented eight-compound mixture with two hard pairs, tuned so about 8 % of the Optimisation 1 design space resolves everything. |
| `simHPLC.m` | Retention in solvent strength and temperature (linear-solvent-strength plus van 't Hoff), migration integrated along the column, plate count from a van Deemter shape, slightly tailed peaks, solvent-front blip, drift and noise. Returns the trace and the ground truth. |
| `simPeaksRow.m` | Peak-picks with the repo's `autofindpeaks` at the saved settings, in the app's data layout. |
| `runSimOptimisation.m` | The loop: design, run, score with the repo objective, suggest, repeat. Logs the truth beside every score. |
| `check_saved_scores.m` | Re-runs the five objectives on the saved peak tables and compares with the saved scores. |
| `test_sim.m` | The checks below. |

## What was checked (Octave 8.4, 2026-11-12; not yet run in MATLAB)

- **Physics.** Isocratic retention matches t0(1+k), and a linear gradient matches the closed-form LSS result within 1 %.
- **Method writer.** All three timetable branches match `RunCond.m`'s rows.
- **Peak picking.** `autofindpeaks` recovers simulated peak positions within 0.01 min and FWHM within 15 %.
- **Inputs that do nothing.** Hotter elutes earlier, and Column/OMPhase/pH change nothing, as declared.
- **Saved scores.** The repo's objectives reproduce the paper's saved scores for all 420 runs (Optimisations 1–9) to 4e-15, with two qualifications:
  - Optimisations 1–3 were saved before the three single objectives gained their log transforms. They match the old forms exactly.
  - `CombinedFunction1` differs by 0.013 and 0.002 in Optimisations 4 and 7, because it normalises by the longest run so far.

## Things the simulator shows about the objectives

1. **The resolution formula reads 0.72 of true resolution.** The objectives take the peak table's width column, which is FWHM (`autofindpeaks` computes area as 1.0646 × height × width), and multiply it by 2.355 as "width at base". Base width is 4σ = 1.699 FWHM. So the objective's "Rs = 2" corresponds to a true Rs of about 2.8. The bias errs toward over-separation and longer methods; it is not unsafe.
2. **The solvent front counts as a peak.** The app calls every objective with `voidTime = 0`, so the filter meant to drop peaks near the solvent front removes nothing. In all 420 saved runs the first peak sits at 0.33–0.39 min, consistent with the front. It sets the critical resolution in only one run, so in their data it mostly adds a constant to the peak count.
3. **Co-eluting peaks can look resolved.** When two compounds merge, the picker sees one peak, and the formula reports good resolution among the peaks it can see. In the demo run, two methods with a true critical Rs of 0.09 and 0.00 scored 1.57 and 1.34 on the repo formula. The peak-count term is the only guard, which is why the weights matter.

In the demo (random search, 32 runs), the objective's five best methods all resolved every compound. Its pick (true Rs 2.05, last peak 4.1 min) was not the truth's pick (Rs 1.61, 3.2 min), because the weights favour resolution over time. That is a weighting choice, not an error.

## Limits

- **Peak widths are conservative.** There is no gradient band compression, so real gradient peaks would be narrower.
- **The parameters are invented.** Void volume, dwell volume, the plate model and the sample are all invented.
- **One channel only.** The trace is labelled 210 nm.
- **Some method inputs are ignored.** Column, mobile phase and pH have no effect.
- **The noise is deterministic.** It is portable rather than statistical, so runs repeat exactly; don't use it to test noise statistics.

The simulator is for learning the optimiser and stress-testing objectives. It is not a model of any real column or sample.

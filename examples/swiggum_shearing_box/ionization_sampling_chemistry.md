# Sampling and hydro-step chemistry experiment

For a **new task-based RHD run**, merge these entries into its existing blocks:

```yaml
TaskBasedRadiationHydrodynamicsSimulation:
  frequency sampling uniform fraction: 0.5
  chemistry every hydro step: true
  chemistry diagnostics: true
  time dependent ionization: true
  advect ionization: true
  do explicit temperature calculation: true
  # Retain detailed line cooling, accelerated by the existing lookup table.
  use cooling tables: false
  line cooling lookup table: true

TemperatureCalculator:
  do temperature calculation: false
```

Leave the photon count, radiation interval, cooling floor, source luminosities,
and other physical parameters unchanged for this comparison. The three new
options default to `0.0`, `false`, and `false`, respectively. There is no new
continuous forcing or smoothing of ion fractions. Advection is unchanged.

## What changes

**Frequency sampling.** Each supported discrete source samples the photon-number
distribution `q = (1-a)*p + a*uniform_frequency`, then assigns weight `p/q`.
`a=0` follows the old sampling and random-number sequence exactly; `a=1` is
uniform frequency sampling. `a=0.5` gives many more hard packets while limiting
the weights to at most two, protecting the abundant softer photons from very
large weights. It does not alter QH0, the spectrum's frequency bounds or cross
sections. Existing path-length ionization/heating estimators use these weights;
diffuse reemission retains them. Printed raw absorbed/escaped packet **counts**
are not luminosity fractions when importance sampling is enabled.

Implemented for WMBasic and Pegase3 spectra, including per-star spectra in
MixedDriving, SwiggumFile and HDF5 source distributions. Other discrete spectrum
types fail clearly if importance sampling is requested. Continuous sources keep
their existing sampler. No static-mode parameter is added by this RHD change.

**Chemistry timing.** Photon shoots remain at the configured radiation interval.
Radiation iterations make trial H/He opacity predictions over one current hydro
timestep, always starting from the same physical state. Those trial states are
discarded after the last iteration; metals are not repeatedly advanced in them.
The final photoionization and photoheating rates are retained in the existing
cell estimator storage, normalized to SI.

After each hydro/advection step, chemistry advances H/He and metals **once** over
that hydro timestep using current density and temperature and the cached rates.
Temperature is then made consistent with the unchanged hydro thermal energy,
and existing explicit heating/cooling runs. This is operator splitting, not a
fully coupled chemistry/cooling integrator: radiation rates remain Eulerian and
fixed between shoots, and temperature remains fixed inside each chemistry solve.

The new mode requires time-dependent ionization, explicit heating/cooling, and
`TemperatureCalculator:do temperature calculation: false`. With radiation off,
cached photo-rates are zero and chemistry is collisional/recombining only.

**Edge cases.** The time-dependent H/He RHS now handles exactly zero H0 and He0
without 0/0. Explicit stages are bounded as an element group so an implicit
highest ion cannot become negative simply because individually floored stages
sum above one. Electron-density expressions protect against negative implicit
He++ from roundoff. These repairs do not change atomic rate tables, the GSL
integrator/tolerances/step cap, or the existing restore-to-initial fallback.

## Diagnostics

When enabled (requires hydro-step chemistry), these are written in the working
directory:

- `chemistry_gsl_counts.tsv`: complete event counts per phase, element and GSL
  status; includes the number of attempted solves as a denominator. Routine
  trace floors/roundoff corrections are separate from failed/restored solves
  and element sums exceeding `1+1e-8`.
- `chemistry_gsl_samples.tsv`: at most 32 informative examples per worker per
  flush. Routine output floors do not consume this quota. Rows include density,
  temperature, requested/reached time, accepted steps, initial/raw output ODE
  vectors, initial/current ions, and radiation rates. The time column is the
  end of the reporting interval, not the exact event time.
- Each subsequent Gadget snapshot gets `ChemistryDiagnostics/Restored` and
  `ChemistryDiagnostics/Corrected`: **complete**, original task-subgrid/cell-order
  masks since the previous output or restart, independent of the example quota. Bits 0–5
  mean physical hydro-step HHe/C/N/O/Ne/S; bits 8–13 mean radiation trials.
  `Corrected` includes harmless floors: a nonzero bit alone is not a pathology.

The count file also records summed worker seconds for physical chemistry and
radiation trials. These are summed parallel worker times, **not elapsed wall
time**. Masks use 8 bytes/original cell (128 MiB for 256^3); there are no masks
for radiation copies. Nothing is allocated for these masks when diagnostics are
off. Counts flush at radiation boundaries, snapshots and normal termination.

Restart layout is unchanged. New-mode restarts refresh the radiation rates
before continuing and do not repeat a completed physical chemistry timestep.
The existing restart mechanism reads parameters from the dump: an old dump does
not acquire these options just by editing an external parameter file. This
experiment is intended to start fresh, without inherited glitter.

## Focused checks

```sh
cmake ..
cmake --build . --target CMacIonize testSamplingHydroChemistry -j 8
ctest -R '^testSamplingHydroChemistry$' --output-on-failure
```

The focused test covers weighted flux, enhanced hard sampling, disabled RNG
compatibility, non-committing radiation trials, cached-rate normalization,
an analytic hydrogen recombination solution, finite fully ionized H/He states,
helium normalization, and diagnostic masks beyond the bounded example limit.
It prints a small chemistry cost/failure comparison as well. A faster or cleaner
full simulation is not guaranteed by passing these checks; this is a combined
practical experiment, not a clean attribution of glitter to one mechanism.

### Local validation (4 September 2026)

- Hydrogen-only and 19-ion executables built with OpenMP and assertions.
- Focused sampling/chemistry/diagnostic tests passed in both configurations.
- The existing last-radiation restart-state test passed.
- An 8^3, two-thread, 19-ion RHD smoke run completed three hydro steps and one
  two-iteration photon shoot: exactly 1,536 physical solves per element,
  1,024 H/He trial solves, finite output, and zero restored solves.
- A stop/write/restart smoke check also counted exactly 1,536 physical solves
  per element in total, with the expected one-time radiation refresh.
- In 300,000 WMBasic draws, natural sampling produced 1,714 packets above
  35.12 eV; `a=0.5` produced 71,631 with weighted equivalent 1,730.93.
- A small fixed-state 50-kyr chemistry benchmark using the actual CHIANTI
  rates found 20 shorter advances cost approximately 1.7 times as much at
  8,000 K and 3 times as much at 10^6 K as one long advance. The latter
  comparison is not equal successful work: all ten long hot-cell H/He solves
  restored their input after hitting the cap; none of the shorter ones did.
  This is not a full-run runtime forecast.

Two unrelated existing tests were also attempted: `testDistributedPhotonSource`
does not compile against the existing three-argument SingleStar constructor;
the MPI `testPhotonBuffer` fails its subgrid-index equality assertion. Their
test/packing sources were not changed by this work. Neither is used to claim
validation of multi-node transport, which is outside this change's scope.

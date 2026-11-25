# DircGR_Calibration

Tools to produce and check the calibration inputs of the ePIC hpDIRC
geometric reconstruction (`DircGR` in [EICrecon](https://github.com/eic/EICrecon)):
the photon-path look-up tables (LUTs) and the chromatic correction. The
finished files go to [eic/epic-data](https://github.com/eic/epic-data)
(`dircgr/`) and come with the geometry (epic `compact/calibrations.xml`).

Everything runs inside `eic-shell` (npsim, eicrecon, ROOT, `$DETECTOR_PATH`).
Analysis code is plain C++ ROOT macros.

## Contents

| path | what |
|---|---|
| `scripts/run_single_hadron.sh` | single π / K / p: npsim → eicrecon (DircGR) → QA page |
| `macros/dircgr_qa.C` | one-page QA of a DircGR output file (θC vs p, resolution, npe and residual vs η, PID) |
| `macros/sprfit.C` | fits DircGR's track resolution parameters `spr`, `sigmaExt` from single-hadron samples |
| `scripts/generateLutData.sh` | LUT photon-gun simulation (npsim), input of `lutgen.C` |
| `lut/guncheck.C` | checks the LUT gun positions against the geometry (run by `generateLutData.sh`) |
| `lut/lutgen.C` | builds `lut_bar<B>.root` (photon paths per MCP pixel) from the gun simulation |
| `scripts/generateChrCorrData.sh` | π+ 3 GeV angle scan, npsim + eicrecon calibration dump, input of `chrcorrfit.C` |
| `chrcorr/chrcorrfit.C` | chromatic-correction fits from the DircGR calibration dump |

## Quick check of the reconstruction: single hadrons

```bash
# 2000 pi+, p flat 0.1-8 GeV, theta flat 30-150 deg, DIRC-only geometry
scripts/run_single_hadron.sh -s pi -n 2000 -p 0.1 -P 8 -o out
# kaons and protons
scripts/run_single_hadron.sh -s k -n 2000 -p 0.1 -P 8 -o out
scripts/run_single_hadron.sh -s p -n 2000 -p 0.1 -P 8 -o out
```

Each run writes `out/<tag>.edm4hep.root` (npsim), `out/<tag>.edm4eic.root`
(eicrecon) and `out/<tag>.pdf`, the QA page:

1. reconstructed θC vs p, with the expected π / K / p curves
2. θC − expected θC of the true species, Gaussian core (mean, σ)
3. track θC resolution and mean offset vs p
4. photons used vs track pseudorapidity at the bar
5. θC − expected vs track pseudorapidity
6. PID: fraction of tracks identified as π / K / p vs p

Example pages (8000 events each, p 0.1-8 GeV, θ 30-150°, DIRC-only
geometry): `macros/examples/dircgr_qa_{pi+,kaon+,proton}_p0.1-8.pdf`.

`scripts/run_single_hadron.sh -h` lists all options (angles, geometry,
threads, seed, ...). Until the DircGR files are in a released geometry,
point to them explicitly: `-l <dir with lut_bar*.root> -c <chrcorr.root>`.
The macro also runs on its own:
`root -l -b -q 'macros/dircgr_qa.C("out/pi+_p0.1-8_th30-150.edm4eic.root")'`.

## Track-resolution parameters

DircGR's hypothesis likelihood uses σ_track² = spr²/N + sigmaExt² (N =
photons used). Fit them from single-hadron outputs (π, K, p together):

```bash
root -l -b -q 'macros/sprfit.C("out/pi+_p0.1-8_th30-150.edm4eic.root out/kaon+_p0.1-8_th30-150.edm4eic.root out/proton_p0.1-8_th30-150.edm4eic.root","sprfit.pdf")'
```

It prints the core σ per bin of N, the fit (N ≥ 25; fewer photons degrade
much faster than 1/√N) and the pull widths with the current and fitted
values; put the fitted values into DircGRConfig.

## LUT files

`lut_bar{0..9}.root`, one per bar (DIRCHits cellID `bar`), tree `dircLUT`
with one entry per optical-photon path: `ipixx`, `ipixy` (MCP pixel),
`OPdirx/y/z` (unit direction at the photon start, in the bar), `OPtime`
(time to the MCP, ns), `OPprob` (fraction of the pixel's photons on this
path), plus `npath`, `ipath`, `nOP`, `OPdirRMS`, `OPtimeRMS` (diagnostics).
The grid is stored next to the tree as `TParameter<int>` (`nPixX` = 107,
`nPixY` = 73, `nBar` = 10, `bar`); DircGR reads it from there and refuses
inconsistent files. Pixels are the DIRCHits MCP cells (3.3 mm):
`ipixx = 53 - cellY`, `ipixy = cellX + 36`.

### Making the LUTs, step by step

1. **Simulate the photon gun** (inside eic-shell, with the geometry the LUT
   is for in `$DETECTOR_PATH`):
   ```bash
   scripts/generateLutData.sh -o lutsim -j 8
   ```
   Per bar, 10 chunks of 2M optical photons (3.18 eV) from a point on the
   bar axis, 2 mm inside the prism-side end face (box 0), isotropic into the
   prism hemisphere (θ 105-180°, flat in cos θ). The gun positions are
   checked against the geometry first (`lut/guncheck.C`: inside `bar_vol`,
   centred, 2 mm from the end face, cellID bar = file index); if the bar
   layout changed, give the new positions with `-X -Y -W -Z`. Output
   `lutsim/bar<B>_c<C>.edm4hep.root`, ~620 MB and ~5 min per chunk, so
   ~60 GB and ~8 CPU-hours for all 10 bars. `-b "4 5"` makes only some bars,
   `-h` lists all options.
2. **Find the paths** (compiled with ACLiC, `+`):
   ```bash
   for b in 0 1 2 3 4 5 6 7 8 9; do
     root -l -b -q "lut/lutgen.C+($b,\"lutsim\",\"LUT\")"
   done
   ```
   Writes `LUT/lut_bar<B>.root` with the grid parameters, ~4 min and
   ~1.1 GB per bar. For each MCP pixel, the photons that hit it are grouped
   by friends-of-friends: two photons are linked if their start directions
   differ by less than 2 mrad and their MCP times by less than 0.2 ns, and a
   group is everything joined by a chain of links. From a point source the
   start direction fixes the whole trajectory, so each group is one optical
   path (one sequence of reflections in lens and prism). Groups with ≥ 20
   photons become LUT paths: mean direction, mean time, OPprob = share of the
   pixel's photons. Settings: `lutgen.C+(bar, indir, outdir, link [mrad],
   tlink [ns], minop)`. Each run prints paths/pixel (~7), the fraction of
   photons in paths (~97%) and the median path width (~3.3 mrad); the file
   also has `hNpath`/`hNop` (paths and photons per pixel).
3. **Check** the new LUT with the single-hadron QA (`run_single_hadron.sh
   -l LUT ...`, above). Then make the chromatic correction with it (next
   section), since that depends on the LUT.

## Chromatic correction

`chrcorr.root`: per photon, DircGR corrects θC linearly in dt/L (arrival
time minus expected, per mm of path), separately for photons going directly
to the prism and via the mirror; the parameters are graphs vs track polar
angle (set `0_0` = π+ 3 GeV). The input is a per-photon dump from DircGR
run with no chromatic correction.

### Making the chromatic correction, step by step

1. **Simulate and dump** (inside eic-shell; eicrecon with the DircGR
   `calibFile` parameter):
   ```bash
   scripts/generateChrCorrData.sh -l LUT -o chrsim -j 8
   ```
   π+ at 3 GeV, 500 events at each of 15 polar angles (25-160°), each at the
   azimuth of a bar centre; npsim, then eicrecon with
   `-PDIRC:DircGR:chrcorrFile=` (off) and
   `-PDIRC:DircGR:calibFile=chrsim/dump/S1_pi3GeV_th<θ>.chrcalib.root`.
   `-l` is the LUT directory (default: the files that come with the
   geometry); use the LUT just made when the geometry changed. ~4 min with
   8 jobs.
2. **Fit**:
   ```bash
   root -l -b -q 'chrcorr/chrcorrfit.C("chrsim/dump","chrcorr.root","chrcorrfit.pdf")'
   ```
   A 4th argument, an older chrcorr file, adds a comparison page.
3. **Check**: rerun the single-hadron QA with `-l LUT -c chrcorr.root` on an
   independent sample (another seed); the θC offset vs angle should be
   within ~1 mrad.

## Updating the files in epic-data

The LUTs and the chromatic correction depend on the bar / lens / prism /
MCP geometry and must be regenerated when it changes. Copy the new files
to `epic-data/dircgr/`, open a pull request there, then update the
epic-data commit hash for `dircgr` in epic `compact/calibrations.xml`.

W.J. Llope (Wayne State University)

License: LGPL-3.0-or-later (see `LICENSE`), as EICrecon.

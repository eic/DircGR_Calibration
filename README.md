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
| `lut/lutmeta.C` | writes the grid (`nPixX`, `nPixY`, `nBar`, `bar`) into the LUT files |
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

`scripts/run_single_hadron.sh -h` lists all options (angles, geometry,
threads, seed, ...). Until the DircGR files are in a released geometry,
point to them explicitly: `-l <dir with lut_bar*.root> -c <chrcorr.root>`.
The macro also runs on its own:
`root -l -b -q 'macros/dircgr_qa.C("out/pi+_p0.1-8_th30-150.edm4eic.root")'`.

## LUT files

`lut_bar{0..9}.root`, one per bar (DIRCHits cellID `bar`), tree `dircLUT`
with one entry per optical-photon path: `ipixx`, `ipixy` (MCP pixel),
`OPdirx/y/z` (unit direction at the photon start, in the bar), `OPtime`
(time to the MCP, ns), plus `npath`, `ipath`, `nOP`, `OPprob`, `OPdirRMS`,
`OPtimeRMS` (diagnostics). The grid is stored next to the tree as
`TParameter<int>` (`nPixX` = 107, `nPixY` = 73, `nBar` = 10, `bar`); DircGR
reads it from there and refuses inconsistent files. For files without it:

```bash
root -l -b -q 'lut/lutmeta.C("LUT")'      # LUT/lut_bar0..9.root, 107 x 73, 10 bars
```

Pixels are the DIRCHits MCP cells (3.3 mm): `ipixx = 53 - cellY`,
`ipixy = cellX + 36`.

LUT production (optical-photon gun from 2 mm inside the prism-side end of
each bar, 3.18 eV, isotropic into the prism hemisphere; path finding by
friends-of-friends grouping): to be added (C macros).

## Chromatic correction

1. pi+ 3 GeV at 15 polar angles (25-160 deg), at the bar-centre azimuth,
   500 events each (npsim, DIRC-only geometry).
2. eicrecon with the per-photon calibration dump and **no** chromatic
   correction:
   ```bash
   eicrecon -Pdd4hep:xml_files=$DETECTOR_PATH/epic_dirc_only.xml \
     -PDIRC:DircGR:lutPath=<LUT dir> -PDIRC:DircGR:chrcorrFile= \
     -PDIRC:DircGR:calibFile=dump/<tag>.chrcalib.root \
     -Ppodio:output_file=/dev/null -Ppodio:output_collections=DircGRRadiatorInfo \
     <tag>.edm4hep.root
   ```
3. Fit:
   `root -l -b -q 'chrcorr/chrcorrfit.C("dump","chrcorr.root","chrcorrfit.pdf")'`
   (a 4th argument, an older chrcorr file, adds a comparison page).

## Updating the files in epic-data

The LUTs and the chromatic correction depend on the bar / lens / prism /
MCP geometry and must be regenerated when it changes. Copy the new files
to `epic-data/dircgr/`, open a pull request there, then update the
epic-data commit hash for `dircgr` in epic `compact/calibrations.xml`.

W.J. Llope (Wayne State University)

#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 W.J. Llope
#
# generateChrCorrData.sh: input of chrcorr/chrcorrfit.C. pi+ at 3 GeV at 15
# polar angles (25-160 deg), each at the azimuth of the bar centre of box 0
# (corrected for the bending at 3 GeV), simulated with npsim, then
# reconstructed with eicrecon (DircGR) with NO chromatic correction and the
# per-photon calibration dump switched on (-PDIRC:DircGR:calibFile).
#
# Run inside eic-shell (needs npsim, eicrecon, $DETECTOR_PATH). The LUTs must
# match the geometry: make them first (generateLutData.sh + lut/lutgen.C)
# when the geometry changed.
#
#   scripts/generateChrCorrData.sh -l LUT -j 8
#   root -l -b -q 'chrcorr/chrcorrfit.C("chrsim/dump","chrcorr.root","chrcorrfit.pdf")'
#
# Options (defaults in brackets):
#   -n events per angle                                         [500]
#   -a angles [deg]          ["25 30 40 50 60 70 80 90 100 110 120 130 140 150 160"]
#   -o outdir                                                   [./chrsim]
#   -l LUT directory (lut_bar0..9.root)   [calibrations/dircgr from the geometry]
#   -g compact file        [$DETECTOR_PATH/epic_dirc_only.xml]
#   -j jobs in parallel                                         [1]
#   -r seed base, seed = base + index of the angle (1..15)      [1000]
#   -x skip npsim if the .edm4hep.root file already exists
#
# Output in outdir, TAG = S1_pi3GeV_th<angle>:
#   TAG.edm4hep.root, TAG.npsim.log, TAG.eicrecon.log
#   dump/TAG.chrcalib.root   (the chrcorrfit.C input)

set -euo pipefail

nev=500; angles="25 30 40 50 60 70 80 90 100 110 120 130 140 150 160"; outdir=./chrsim
lut=""; compact=""; npar=1; seedbase=1000; skipsim=0

usage() { sed -n '5,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts "n:a:o:l:g:j:r:xh" opt; do
  case $opt in
    n) nev=$OPTARG ;; a) angles=$OPTARG ;; o) outdir=$OPTARG ;; l) lut=$OPTARG ;;
    g) compact=$OPTARG ;; j) npar=$OPTARG ;; r) seedbase=$OPTARG ;; x) skipsim=1 ;;
    h) usage 0 ;; *) usage 1 ;;
  esac
done

for cmd in npsim eicrecon; do
  command -v $cmd >/dev/null || { echo "$cmd not found: run inside eic-shell"; exit 1; }
done
: "${DETECTOR_PATH:?DETECTOR_PATH not set: run inside eic-shell}"
compact=${compact:-$DETECTOR_PATH/epic_dirc_only.xml}
[ -f "$compact" ] || { echo "compact file $compact not found"; exit 1; }
if [ -n "$lut" ]; then
  lut=$(cd "$lut" && pwd)
  for b in 0 1 2 3 4 5 6 7 8 9; do
    [ -f "$lut/lut_bar$b.root" ] || { echo "no $lut/lut_bar$b.root"; exit 1; }
  done
fi
mkdir -p "$outdir/dump"
outdir=$(cd "$outdir" && pwd)

# azimuth of the bar centre (box 0, bar 5 side) for a 3 GeV pi+ from the
# origin, vs polar angle: table, linear interpolation
phi_of_theta() {
  awk -v t="$1" 'BEGIN{
    n = split("22.65 30 50 70 86 88 90 92 94 110 125 130 150 164.3", T, " ");
    split("350.76 351.95 353.93 354.73 354.941 354.945 354.95 354.945 354.941 354.70 354.15 353.84 351.68 348.85", P, " ");
    if (t <= T[1]) { printf "%.3f", P[1]; exit }
    if (t >= T[n]) { printf "%.3f", P[n]; exit }
    for (i = 1; i < n; i++) if (t >= T[i] && t <= T[i+1]) {
      printf "%.3f", P[i] + (t - T[i]) * (P[i+1] - P[i]) / (T[i+1] - T[i]); exit }
  }'
}

run_angle() {
  local th=$1 seed=$2
  local tag=S1_pi3GeV_th$th ph
  local sim=$outdir/$tag.edm4hep.root dump=$outdir/dump/$tag.chrcalib.root
  ph=$(phi_of_theta "$th")
  if [ $skipsim = 1 ] && [ -f "$sim" ]; then
    echo "-- $tag: using existing $sim"
  else
    echo "-- $tag: npsim $nev pi+ 3 GeV, theta $th, phi $ph deg, seed $seed"
    npsim --runType batch --printLevel WARNING --random.seed "$seed" \
      --compactFile "$compact" -N "$nev" --enableGun --gun.particle pi+ \
      --gun.momentumMin "3*GeV" --gun.momentumMax "3*GeV" \
      --gun.thetaMin "${th}*deg" --gun.thetaMax "${th}*deg" \
      --gun.phiMin "${ph}*deg" --gun.phiMax "${ph}*deg" --gun.distribution uniform \
      --gun.position "0*cm,0*cm,0*cm" --outputFile "$sim" > "$outdir/$tag.npsim.log" 2>&1 \
      || { echo "npsim failed: $outdir/$tag.npsim.log"; return 1; }
  fi
  echo "-- $tag: eicrecon -> $dump"
  local opts=(-Pdd4hep:xml_files="$compact" -Pnthreads=1
              -PDIRC:DircGR:chrcorrFile= -PDIRC:DircGR:calibFile="$dump"
              -Ppodio:output_file=/dev/null -Ppodio:output_collections=DircGRRadiatorInfo)
  [ -n "$lut" ] && opts+=(-PDIRC:DircGR:lutPath="$lut")
  rm -f "$dump"
  (cd "$outdir" && eicrecon "${opts[@]}" "$sim") > "$outdir/$tag.eicrecon.log" 2>&1 \
    || { echo "eicrecon failed: $outdir/$tag.eicrecon.log"; return 1; }
  [ -s "$dump" ] || { echo "no $dump: does this eicrecon have DIRC:DircGR:calibFile?"; return 1; }
}

echo "== chromatic-correction input: pi+ 3 GeV, $nev events at [$angles] deg -> $outdir"
echo "   geometry $compact, LUT ${lut:-calibrations/dircgr (geometry default)}"
pids=(); i=0
for th in $angles; do
  i=$((i + 1))
  run_angle "$th" $((seedbase + i)) &
  pids+=($!)
  while [ "$(jobs -rp | wc -l)" -ge "$npar" ]; do sleep 2; done
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ $fail = 0 ] || { echo "== some angles FAILED"; exit 1; }
echo "== done. Next: root -l -b -q '$(cd "$(dirname "$0")/.." && pwd)/chrcorr/chrcorrfit.C(\"$outdir/dump\",\"chrcorr.root\",\"chrcorrfit.pdf\")'"

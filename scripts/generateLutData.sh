#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 W.J. Llope
#
# generateLutData.sh: LUT photon-gun simulation (npsim), the input of
# lut/lutgen.C. Per bar: optical photons at 3.18 eV from a point at the bar
# centre, 2 mm inside the prism-side end face of the bar (box 0), isotropic
# into the prism hemisphere (theta 105-180 deg, flat in cos(theta)).
# The gun positions are checked against the geometry first (lut/guncheck.C).
#
# Run inside eic-shell (needs npsim, root, $DETECTOR_PATH).
#
#   scripts/generateLutData.sh                      # all 10 bars, 10 x 2M photons each
#   scripts/generateLutData.sh -b "4 5" -c 2 -j 2   # bars 4, 5, 2 chunks, 2 at a time
#
# Options (defaults in brackets):
#   -b bars                list of cellID bars                 ["0 1 2 3 4 5 6 7 8 9"]
#   -c chunks per bar                                           [10]
#   -n photons per chunk                                        [2000000]
#   -o outdir                                                   [./lutsim]
#   -g compact file        [$DETECTOR_PATH/epic_dirc_only.xml]
#   -j npsim jobs in parallel                                   [1]
#   -r seed base, seed = base + 100*bar + chunk                 [100000]
#   -x skip chunks whose output already exists
#   -k skip the gun-position check (studies with a deliberately different gun)
#   gun position [mm], bar b at (X, Y0 - W*b, Z):
#   -X x  -Y y0  -W pitch  -Z z                [765.875, 158.175, 35.15, -2725.55]
#
# Output: outdir/bar<B>_c<C>.edm4hep.root (~620 MB per 2M photons, ~5 min),
# logs outdir/bar<B>_c<C>.log. Then, per bar:
#   root -l -b -q 'lut/lutgen.C+(B,"lutsim","LUT")'

set -euo pipefail

bars="0 1 2 3 4 5 6 7 8 9"; nchunk=10; nphot=2000000; outdir=./lutsim; compact=""
npar=1; seedbase=100000; skip=0; check=1
X=765.875; Y0=158.175; W=35.15; Z=-2725.55

usage() { sed -n '5,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts "b:c:n:o:g:j:r:xkX:Y:W:Z:h" opt; do
  case $opt in
    b) bars=$OPTARG ;; c) nchunk=$OPTARG ;; n) nphot=$OPTARG ;; o) outdir=$OPTARG ;;
    g) compact=$OPTARG ;; j) npar=$OPTARG ;; r) seedbase=$OPTARG ;; x) skip=1 ;; k) check=0 ;;
    X) X=$OPTARG ;; Y) Y0=$OPTARG ;; W) W=$OPTARG ;; Z) Z=$OPTARG ;;
    h) usage 0 ;; *) usage 1 ;;
  esac
done

for cmd in npsim root; do
  command -v $cmd >/dev/null || { echo "$cmd not found: run inside eic-shell"; exit 1; }
done
: "${DETECTOR_PATH:?DETECTOR_PATH not set: run inside eic-shell}"
compact=${compact:-$DETECTOR_PATH/epic_dirc_only.xml}
[ -f "$compact" ] || { echo "compact file $compact not found"; exit 1; }
here=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)

echo "== LUT gun: bars [$bars], $nchunk x $nphot photons each -> $outdir"
echo "   geometry $compact"
if [ $check = 1 ]; then
  root -l -b -q "$here/lut/guncheck.C(\"$compact\",$X,$Y0,$W,$Z)" > "$outdir/guncheck.log" 2>&1 \
    || { grep -E "^bar|guncheck" "$outdir/guncheck.log"; echo "gun positions do not match the geometry (-X -Y -W -Z)"; exit 1; }
  grep -E "guncheck" "$outdir/guncheck.log"
else
  echo "   gun-position check SKIPPED (-k)"
fi

run_chunk() {
  local b=$1 c=$2
  local y seed out
  y=$(awk -v y0="$Y0" -v w="$W" -v b="$b" 'BEGIN{printf "%.3f", y0 - w*b}')
  seed=$((seedbase + 100 * b + c))
  out=$outdir/bar${b}_c${c}.edm4hep.root
  if [ $skip = 1 ] && [ -f "$out" ]; then echo "-- bar $b chunk $c: exists"; return 0; fi
  echo "-- bar $b chunk $c: gun ($X, $y, $Z) mm, seed $seed"
  npsim --runType batch --compactFile "$compact" -G -N "$nphot" \
    --enableGun --printLevel WARNING --random.seed "$seed" \
    --gun.multiplicity 1 --gun.particle opticalphoton --gun.energy "3.18*eV" \
    --gun.phiMin "0*deg" --gun.phiMax "360*deg" --gun.thetaMin "105*deg" \
    --gun.thetaMax "180*deg" --gun.distribution "cos(theta)" \
    --gun.position "${X}*mm,${y}*mm,${Z}*mm" \
    --part.userParticleHandler='' --outputFile "$out.tmp.root" > "$outdir/bar${b}_c${c}.log" 2>&1 \
    || { echo "npsim failed: $outdir/bar${b}_c${c}.log"; rm -f "$out.tmp.root"; return 1; }
  mv "$out.tmp.root" "$out"
}

pids=()
for b in $bars; do
  for ((c = 0; c < nchunk; c++)); do
    run_chunk "$b" "$c" &
    pids+=($!)
    while [ "$(jobs -rp | wc -l)" -ge "$npar" ]; do sleep 2; done
  done
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ $fail = 0 ] || { echo "== some chunks FAILED"; exit 1; }
echo "== done. Next: for b in $bars; do root -l -b -q \"$here/lut/lutgen.C+(\$b,\\\"$outdir\\\",\\\"LUT\\\")\"; done"

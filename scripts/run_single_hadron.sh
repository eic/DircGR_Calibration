#!/bin/bash
# run_single_hadron.sh: single charged hadrons through the ePIC DIRC
# (npsim) and the EICrecon DircGR reconstruction, then the one-page QA plot
# (macros/dircgr_qa.C).
#
# Run inside eic-shell (needs npsim, eicrecon, root, $DETECTOR_PATH).
#
#   scripts/run_single_hadron.sh -s pi+ -n 2000 -p 0.1 -P 8
#   scripts/run_single_hadron.sh -s kaon- -n 500 -t 60 -T 90 -o out_k
#
# Options (defaults in brackets):
#   -s species   pi+ pi- kaon+ kaon- proton anti_proton, or pi / k / p
#                (= pi+ / kaon+ / proton)                       [pi+]
#   -n events                                                   [1000]
#   -p pmin  -P pmax       momentum, flat [GeV]                 [0.1, 8]
#   -t thmin -T thmax      polar angle, flat in theta [deg]     [30, 150]
#   -f phimin -F phimax    azimuth [deg]                        [0, 360]
#   -o outdir                                                   [./out]
#   -g compact file        [$DETECTOR_PATH/epic_dirc_only.xml]
#                          (DIRC only, fast; $DETECTOR_PATH/epic.xml = full ePIC)
#   -j eicrecon threads                                         [1]
#   -r seed                                                     [1]
#   -l LUT directory (lut_bar0..9.root)  [calibrations/dircgr from the geometry]
#   -c chromatic-correction file         [calibrations/dircgr/chrcorr.root]
#   -x skip npsim if the .edm4hep.root file already exists
#   -q no QA plot
#
# Output in outdir, TAG = <species>_p<pmin>-<pmax>_th<thmin>-<thmax>:
#   TAG.edm4hep.root   npsim
#   TAG.edm4eic.root   eicrecon (MCParticles + DIRC/DircGR collections)
#   TAG.pdf            QA page (thetaC vs p, residuals, npe, PID)
#   TAG.npsim.log, TAG.eicrecon.log

set -euo pipefail

species=pi+; nev=1000; pmin=0.1; pmax=8; thmin=30; thmax=150; phimin=0; phimax=360
outdir=./out; compact=""; nthr=1; seed=1; lut=""; chr=""; skipsim=0; qa=1

usage() { sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while getopts "s:n:p:P:t:T:f:F:o:g:j:r:l:c:xqh" opt; do
  case $opt in
    s) species=$OPTARG ;; n) nev=$OPTARG ;;
    p) pmin=$OPTARG ;; P) pmax=$OPTARG ;;
    t) thmin=$OPTARG ;; T) thmax=$OPTARG ;;
    f) phimin=$OPTARG ;; F) phimax=$OPTARG ;;
    o) outdir=$OPTARG ;; g) compact=$OPTARG ;; j) nthr=$OPTARG ;; r) seed=$OPTARG ;;
    l) lut=$OPTARG ;; c) chr=$OPTARG ;; x) skipsim=1 ;; q) qa=0 ;;
    h) usage 0 ;; *) usage 1 ;;
  esac
done

case $species in
  pi) species=pi+ ;; k|K|kaon) species=kaon+ ;; p|proton) species=proton ;;
  pi+|pi-|kaon+|kaon-|anti_proton) ;;
  *) echo "unknown species '$species' (pi+ pi- kaon+ kaon- proton anti_proton)"; exit 1 ;;
esac

for cmd in npsim eicrecon root; do
  command -v $cmd >/dev/null || { echo "$cmd not found: run inside eic-shell"; exit 1; }
done
: "${DETECTOR_PATH:?DETECTOR_PATH not set: run inside eic-shell}"
compact=${compact:-$DETECTOR_PATH/epic_dirc_only.xml}
[ -f "$compact" ] || { echo "compact file $compact not found"; exit 1; }

here=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)
tag=${species}_p${pmin}-${pmax}_th${thmin}-${thmax}
sim=$outdir/$tag.edm4hep.root
rec=$outdir/$tag.edm4eic.root

echo "== $tag: $nev x $species, p $pmin-$pmax GeV, theta $thmin-$thmax deg, phi $phimin-$phimax deg"
echo "   geometry $compact"

# 1. npsim. The gun sits at the origin; the DIRC bars are at r ~ 77 cm.
if [ $skipsim = 1 ] && [ -f "$sim" ]; then
  echo "-- npsim: using existing $sim"
else
  echo "-- npsim -> $sim"
  npsim --runType batch --printLevel WARNING --random.seed "$seed" \
    --compactFile "$compact" -N "$nev" \
    --enableGun --gun.particle "$species" \
    --gun.momentumMin "${pmin}*GeV" --gun.momentumMax "${pmax}*GeV" \
    --gun.thetaMin "${thmin}*deg" --gun.thetaMax "${thmax}*deg" \
    --gun.phiMin "${phimin}*deg" --gun.phiMax "${phimax}*deg" \
    --gun.distribution uniform --gun.position "0*cm,0*cm,0*cm" \
    --outputFile "$sim" > "$outdir/$tag.npsim.log" 2>&1 \
    || { echo "npsim failed, see $outdir/$tag.npsim.log"; tail -20 "$outdir/$tag.npsim.log"; exit 1; }
fi

# 2. eicrecon: DIRC digitization (DircGRdig) + truth bar incidence + DircGR.
# The calibration files come with the geometry (calibrations/dircgr, fetched
# into the run directory); -l / -c override them.
echo "-- eicrecon -> $rec"
opts=(-Pdd4hep:xml_files="$compact" -Pnthreads="$nthr"
      -Ppodio:output_file="$rec"
      -Ppodio:output_collections=MCParticles,DIRCRawHits,DircTruthTrackSegments,DircGRParticleIDs,DircGRRadiatorInfo)
[ -n "$lut" ] && opts+=(-PDIRC:DircGR:lutPath="$lut")
[ -n "$chr" ] && opts+=(-PDIRC:DircGR:chrcorrFile="$chr")
(cd "$outdir" && eicrecon "${opts[@]}" "$sim") > "$outdir/$tag.eicrecon.log" 2>&1 \
  || { echo "eicrecon failed, see $outdir/$tag.eicrecon.log"; tail -20 "$outdir/$tag.eicrecon.log"; exit 1; }
if grep -q "no LUT\|no chromatic correction" "$outdir/$tag.eicrecon.log"; then
  echo "WARNING: DircGR ran without LUT or chromatic correction:"
  grep "no LUT\|no chromatic correction" "$outdir/$tag.eicrecon.log" | sort -u | sed 's/^/   /'
  echo "   (use -l <LUT dir> -c <chrcorr.root> until the files come with the geometry)"
fi

# 3. QA page
if [ $qa = 1 ]; then
  echo "-- QA -> $outdir/$tag.pdf"
  root -l -b -q "$here/macros/dircgr_qa.C(\"$rec\",\"$outdir/$tag.pdf\")" \
    > "$outdir/$tag.qa.log" 2>&1 || { echo "QA macro failed, see $outdir/$tag.qa.log"; exit 1; }
  grep -E "^(tracks|thetaC|PID)" "$outdir/$tag.qa.log" || true
fi
echo "== done: $outdir/$tag.*"

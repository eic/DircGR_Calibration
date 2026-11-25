// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 W.J. Llope
//
// guncheck.C -- check the LUT photon-gun positions against the geometry.
// For each bar b the gun sits at (x, y0 - pitch*b, z) [mm]. It must be
//   - inside a volume named bar_vol, centred across the bar (in x and y),
//   - `depth` mm from the prism-side (-z) end face of the bar,
//   - in DIRCHits cellID module 0 and cellID bar b (the LUT file index).
// Prints one line per bar and returns the number of bad bars.
//
//   root -l -b -q 'lut/guncheck.C("'$DETECTOR_PATH'/epic_dirc_only.xml")'
//   root -l -b -q 'lut/guncheck.C("epic_dirc_only.xml", 765.875, 158.175, 35.15, -2725.55, 2.0)'

#include <DD4hep/Detector.h>
#include <DD4hep/Readout.h>
#include <DDRec/CellIDPositionConverter.h>
#include <TGeoManager.h>
#include <TGeoNavigator.h>
#include <TString.h>

#include <cmath>
#include <cstdio>

R__LOAD_LIBRARY(libDDCore)
R__LOAD_LIBRARY(libDDRec)

int guncheck(const char* compact, double x = 765.875, double y0 = 158.175, double pitch = 35.15,
             double z = -2725.55, double depth = 2.0, int nBar = 10) {
  auto& det = dd4hep::Detector::getInstance();
  det.fromXML(compact);
  TGeoManager* g = gGeoManager;
  dd4hep::rec::CellIDPositionConverter conv(det);
  auto* dec = det.readout("DIRCHits").idSpec().decoder();
  int nbad  = 0;
  // geometry in cm (TGeo / dd4hep units), printout in mm
  for (int b = 0; b < nBar; b++) {
    const double y = y0 - pitch * b;
    g->SetCurrentPoint(x / 10, y / 10, z / 10);
    TGeoNode* n          = g->FindNode();
    const TString name   = n ? n->GetVolume()->GetName() : "none";
    auto dist            = [&](double ux, double uy, double uz) {
      g->SetCurrentPoint(x / 10, y / 10, z / 10);
      g->SetCurrentDirection(ux, uy, uz);
      g->FindNode();
      g->FindNextBoundary(1e4);
      return g->GetStep() * 10;
    };
    const double dx = (dist(1, 0, 0) - dist(-1, 0, 0)) / 2;
    const double dy = (dist(0, 1, 0) - dist(0, -1, 0)) / 2;
    const double dz = dist(0, 0, -1);
    long long mod = -1, cbar = -1;
    const auto cid = conv.cellID(dd4hep::Position(x / 10, y / 10, z / 10));
    if (cid != 0) {
      mod  = dec->get(cid, "module");
      cbar = dec->get(cid, "bar");
    }
    const bool ok = name == "bar_vol" && std::fabs(dx) < 1e-3 && std::fabs(dy) < 1e-3 &&
                    std::fabs(dz - depth) < 1e-3 && mod == 0 && cbar == b;
    if (!ok) nbad++;
    printf("bar %d gun (%.3f, %.3f, %.3f) mm: %-8s off-centre x %+.4f y %+.4f, %.4f mm from end "
           "face, cellID module %lld bar %lld  %s\n",
           b, x, y, z, name.Data(), dx, dy, dz, mod, cbar, ok ? "OK" : "BAD");
  }
  printf(nbad ? "guncheck: %d BAD bar(s)\n" : "guncheck: all OK\n", nbad);
  return nbad;
}

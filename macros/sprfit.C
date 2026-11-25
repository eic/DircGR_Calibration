// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 W.J. Llope
//
// sprfit.C: fit the DircGR track-resolution parameters spr and sigmaExt
//
//   sigma_track^2 = spr^2 / N + sigmaExt^2      (mrad, N = photons used)
//
// from single-hadron samples (scripts/run_single_hadron.sh). Per file, the
// generated particle is matched to its bar-entry segment (closest |p|,
// within 5%), dtheta = thetaC - expected(true species, n) is taken, and the
// file's Gaussian-core mean is subtracted (species-dependent offsets are a
// calibration issue, not resolution). In bins of N the Gaussian core sigma
// of dtheta is fitted, then sigma(N) with the formula above. The fitted
// values go into DircGRConfig (spr, sigmaExt).
//
//   root -l -b -q 'sprfit.C("out/pi+_p0.1-8_th30-150.edm4eic.root out/kaon+_p0.1-8_th30-150.edm4eic.root","sprfit.pdf")'
//
// The fit uses N >= nmin (default 25): fewer photons degrade much faster
// than 1/sqrt(N) (peak finder), see the printed per-bin table.
// Also prints the pull width (dtheta / sigma_track with the CURRENT
// DircGR config values spr0, sigmaExt0) for comparison.

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <vector>

namespace sfit {

double thetaExp(double p, double m, double n) {
  const double beta = p / std::sqrt(p * p + m * m);
  return (n * beta > 1) ? 1000. * std::acos(1. / (n * beta)) : 0.;
}

double massOf(int pdg) {
  switch (std::abs(pdg)) {
  case 211: return 0.13957039;
  case 321: return 0.493677;
  case 2212: return 0.93827208;
  default: return -1;
  }
}

std::vector<double> leaf(TTree* t, const char* name) {
  std::vector<double> v;
  TLeaf* l = t->GetLeaf(name);
  if (l != nullptr) {
    for (int k = 0; k < l->GetLen(); k++) {
      v.push_back(l->GetValue(k));
    }
  }
  return v;
}

// Gaussian core: fit, then refit twice within +-2 sigma; nullptr if too few
TF1* fitCore(TH1* h, int minEntries = 30) {
  if (h->GetEntries() < minEntries) {
    return nullptr;
  }
  TF1* f = new TF1(TString::Format("f_%s", h->GetName()), "gaus", h->GetXaxis()->GetXmin(),
                   h->GetXaxis()->GetXmax());
  f->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
  h->Fit(f, "QN0");
  for (int it = 0; it < 2; it++) {
    const double m = f->GetParameter(1), s = std::abs(f->GetParameter(2));
    f->SetRange(m - 2 * s, m + 2 * s);
    h->Fit(f, "QRN0");
  }
  return f;
}

struct Trk {
  double dth; // thetaC - expected [mrad]
  double npe; // photons used
};

// matched primary tracks with a thetaC and an expected angle
std::vector<Trk> readFile(const char* fname, double n) {
  std::vector<Trk> out;
  TFile f(fname);
  auto* t = f.Get<TTree>("events");
  if (t == nullptr) {
    printf("no events tree in %s\n", fname);
    return out;
  }
  t->SetBranchStatus("*", false);
  for (const char* b : {"MCParticles*", "DircTruthTrackSegments*", "_DircTruthTrackSegments_points*",
                        "DircGRRadiatorInfo*"}) {
    t->SetBranchStatus(b, true);
  }
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    const auto st  = leaf(t, "MCParticles.generatorStatus");
    const auto pdg = leaf(t, "MCParticles.PDG");
    const auto px = leaf(t, "MCParticles.momentum.x"), py = leaf(t, "MCParticles.momentum.y"),
               pz = leaf(t, "MCParticles.momentum.z");
    double pgen = -1, m = -1;
    for (size_t k = 0; k < st.size(); k++) {
      if (st[k] == 1) {
        pgen = std::sqrt(px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k]);
        m    = massOf(static_cast<int>(pdg[k]));
        break;
      }
    }
    if (pgen <= 0 || m < 0) {
      continue;
    }
    const auto pb  = leaf(t, "DircTruthTrackSegments.points_begin");
    const auto pe  = leaf(t, "DircTruthTrackSegments.points_end");
    const auto sx  = leaf(t, "_DircTruthTrackSegments_points.momentum.x");
    const auto sy  = leaf(t, "_DircTruthTrackSegments_points.momentum.y");
    const auto sz  = leaf(t, "_DircTruthTrackSegments_points.momentum.z");
    const auto ang = leaf(t, "DircGRRadiatorInfo.angle");
    const auto npe = leaf(t, "DircGRRadiatorInfo.npe");
    int best       = -1;
    double pbest = 0, dbest = 1e9;
    for (size_t s = 0; s < pb.size(); s++) {
      if (pe[s] <= pb[s]) {
        continue;
      }
      const auto k   = static_cast<size_t>(pb[s]);
      const double p = std::sqrt(sx[k] * sx[k] + sy[k] * sy[k] + sz[k] * sz[k]);
      if (std::abs(p - pgen) < dbest) {
        dbest = std::abs(p - pgen);
        best  = static_cast<int>(s);
        pbest = p;
      }
    }
    if (best < 0 || dbest > 0.05 * pgen || static_cast<size_t>(best) >= ang.size()) {
      continue;
    }
    const double thE = thetaExp(pbest, m, n);
    if (ang[best] <= 0 || thE <= 0) {
      continue;
    }
    out.push_back({ang[best] - thE, npe[best]});
  }
  return out;
}

} // namespace sfit

void sprfit(const char* files, const char* pdf = "sprfit.pdf", double n = 1.4738,
            double spr0 = 6.1, double sigmaExt0 = 0.86, double nmin = 25) {
  using namespace sfit;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // N bins (photons used)
  const std::vector<double> nb = {3, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300};
  const int nbin               = static_cast<int>(nb.size()) - 1;
  auto* h2   = new TH2D("h2", "", nbin, nb.data(), 200, -10, 10);
  auto* hPu0 = new TH1D("hPu0", "pull, current spr/sigmaExt;#Delta#theta / #sigma_{track};tracks", 120, -6, 6);
  auto* hPu1 = new TH1D("hPu1", "pull, fitted spr/sigmaExt;#Delta#theta / #sigma_{track};tracks", 120, -6, 6);

  std::vector<Trk> all;
  TObjArray* names = TString(files).Tokenize(" ,");
  for (int i = 0; i < names->GetEntries(); i++) {
    const char* fn = static_cast<TObjString*>(names->At(i))->GetString().Data();
    auto trk       = readFile(fn, n);
    // subtract this file's (species') Gaussian-core mean
    TH1D hm("hm", "", 160, -20, 20);
    for (const auto& t : trk) {
      hm.Fill(t.dth);
    }
    TF1* fm         = fitCore(&hm);
    const double mu = fm ? fm->GetParameter(1) : 0.;
    printf("%s: %zu tracks, core mean %.2f mrad subtracted\n", fn, trk.size(), mu);
    for (auto& t : trk) {
      t.dth -= mu;
      all.push_back(t);
      h2->Fill(t.npe, t.dth);
    }
  }

  // sigma(N)
  auto* g = new TGraphErrors();
  for (int b = 1; b <= nbin; b++) {
    TH1D* hy = h2->ProjectionY(TString::Format("hy%d", b), b, b);
    TF1* fs  = fitCore(hy);
    if (fs == nullptr) {
      continue;
    }
    // mean N in the bin
    double sn = 0, cnt = 0;
    for (const auto& t : all) {
      if (t.npe >= nb[b - 1] && t.npe < nb[b]) {
        sn += t.npe;
        cnt++;
      }
    }
    printf("N %5.0f-%5.0f: %6.0f tracks, <N> %6.1f, core sigma %.2f +- %.2f mrad\n", nb[b - 1], nb[b],
           cnt, sn / cnt, std::abs(fs->GetParameter(2)), fs->GetParError(2));
    g->AddPoint(sn / cnt, std::abs(fs->GetParameter(2)));
    g->SetPointError(g->GetN() - 1, 0, fs->GetParError(2));
  }
  // fit range N >= nmin: below ~25 photons the peak finder degrades much
  // faster than 1/sqrt(N) and the formula does not apply
  TF1* fres = new TF1("fres", "sqrt([0]*[0]/x + [1]*[1])", nmin, 300);
  fres->SetParameters(spr0, sigmaExt0);
  fres->SetParNames("spr", "sigmaExt");
  g->Fit(fres, "QR");
  const double spr = std::abs(fres->GetParameter(0)), ext = std::abs(fres->GetParameter(1));

  // pulls with the old and the new parameters (N >= nmin)
  for (const auto& t : all) {
    if (t.npe < nmin) {
      continue;
    }
    hPu0->Fill(t.dth / std::sqrt(spr0 * spr0 / t.npe + sigmaExt0 * sigmaExt0));
    hPu1->Fill(t.dth / std::sqrt(spr * spr / t.npe + ext * ext));
  }
  TF1* fp0 = fitCore(hPu0);
  TF1* fp1 = fitCore(hPu1);

  auto* c = new TCanvas("c", "sprfit", 1500, 500);
  c->Divide(3, 1);
  c->cd(1);
  gPad->SetLogx();
  g->SetTitle("track #theta_{C} resolution vs photons used;N photons used;#sigma(#Delta#theta_{C}) core [mrad]");
  g->SetMarkerStyle(20);
  g->Draw("ap");
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.045);
  tx.DrawLatex(0.4, 0.82, TString::Format("spr = %.2f #pm %.2f mrad", spr, fres->GetParError(0)));
  tx.DrawLatex(0.4, 0.76, TString::Format("#sigma_{ext} = %.2f #pm %.2f mrad", ext, fres->GetParError(1)));
  c->cd(2);
  hPu0->Draw();
  if (fp0) {
    tx.DrawLatex(0.15, 0.82, TString::Format("#sigma = %.2f (spr %.2f, #sigma_{ext} %.2f)", std::abs(fp0->GetParameter(2)), spr0, sigmaExt0));
  }
  c->cd(3);
  hPu1->Draw();
  if (fp1) {
    tx.DrawLatex(0.15, 0.82, TString::Format("#sigma = %.2f (spr %.2f, #sigma_{ext} %.2f)", std::abs(fp1->GetParameter(2)), spr, ext));
  }
  c->SaveAs(pdf);

  printf("fit: spr = %.2f +- %.2f mrad, sigmaExt = %.2f +- %.2f mrad (%zu tracks)\n", spr,
         fres->GetParError(0), ext, fres->GetParError(1), all.size());
  printf("pull core sigma: current (%.2f, %.2f) %.2f, fitted %.2f\n", spr0, sigmaExt0,
         fp0 ? std::abs(fp0->GetParameter(2)) : 0., fp1 ? std::abs(fp1->GetParameter(2)) : 0.);
}

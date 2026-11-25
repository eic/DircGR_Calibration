// dircgr_qa.C: one-page QA of the EICrecon DircGR reconstruction for
// single-hadron samples (scripts/run_single_hadron.sh).
//
// Input: eicrecon output with MCParticles, DircTruthTrackSegments,
// DircGRRadiatorInfo and DircGRParticleIDs. Per event, the generated
// particle (generatorStatus 1) is matched to the bar-entry segment with the
// closest |p| (within 5%); delta rays and other secondaries are ignored.
//
// Page (3 x 2):
//   1. thetaC vs p, with the expected pi / K / p curves
//   2. thetaC - expected (true species), Gaussian core fit
//   3. track thetaC resolution and mean offset vs p
//   4. photons used (npe) vs track polar angle at the bar
//   5. thetaC - expected vs track polar angle
//   6. PID: fraction identified as pi / K / p vs p (best hypothesis)
//
//   root -l -b -q 'dircgr_qa.C("out/pi+_p1-6_th30-150.edm4eic.root")'
//   root -l -b -q 'dircgr_qa.C("in.edm4eic.root","qa.pdf",1.4738)'
//
// n = refractive index used for the expected angles (DircGR's
// refractiveIndex, default 1.4738).

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace dqa {

const double kMass[3] = {0.13957039, 0.493677, 0.93827208}; // pi, K, p [GeV]
const int kPdg[3]     = {211, 321, 2212};
const char* kName[3]  = {"#pi", "K", "p"};
const char* kAscii[3]  = {"pi", "K", "p"}; // for the printout
const int kColor[3]   = {kBlue + 1, kGreen + 2, kRed + 1};

// expected Cherenkov angle [mrad]; 0 below threshold
double thetaExp(double p, double m, double n) {
  const double beta = p / std::sqrt(p * p + m * m);
  return (n * beta > 1) ? 1000. * std::acos(1. / (n * beta)) : 0.;
}

int species(int pdg) {
  for (int k = 0; k < 3; k++) {
    if (std::abs(pdg) == kPdg[k]) {
      return k;
    }
  }
  return -1;
}

// all values of one leaf in the current entry (float or double members)
std::vector<double> leaf(TTree* t, const char* name) {
  std::vector<double> v;
  TLeaf* l = t->GetLeaf(name);
  if (l == nullptr) {
    return v;
  }
  for (int k = 0; k < l->GetLen(); k++) {
    v.push_back(l->GetValue(k));
  }
  return v;
}

// Gaussian core: fit, then refit within +-2 sigma
TF1* fitCore(TH1* h) {
  if (h->GetEntries() < 20) {
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

} // namespace dqa

void dircgr_qa(const char* infile, const char* pdf = "", double n = 1.4738) {
  using namespace dqa;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetTitleSize(0.05, "t");

  TFile f(infile);
  auto* t = f.Get<TTree>("events");
  if (t == nullptr) {
    printf("no events tree in %s\n", infile);
    return;
  }
  t->SetBranchStatus("*", false);
  for (const char* b : {"MCParticles*", "DircTruthTrackSegments*", "_DircTruthTrackSegments_points*",
                        "DircGRRadiatorInfo*", "DircGRParticleIDs*", "_DircGRParticleIDs_hypotheses*"}) {
    t->SetBranchStatus(b, true);
  }

  // first pass: momentum range and true species
  double pmin = 1e9, pmax = 0;
  int truth = -1;
  const Long64_t nent = t->GetEntries();
  for (Long64_t i = 0; i < nent; i++) {
    t->GetEntry(i);
    const auto st  = leaf(t, "MCParticles.generatorStatus");
    const auto pdg = leaf(t, "MCParticles.PDG");
    const auto px = leaf(t, "MCParticles.momentum.x"), py = leaf(t, "MCParticles.momentum.y"),
               pz = leaf(t, "MCParticles.momentum.z");
    for (size_t k = 0; k < st.size(); k++) {
      if (st[k] == 1) {
        const double p = std::sqrt(px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k]);
        pmin           = std::min(pmin, p);
        pmax           = std::max(pmax, p);
        if (truth < 0) {
          truth = species(static_cast<int>(pdg[k]));
        }
        break;
      }
    }
  }
  if (truth < 0 || pmax <= 0) {
    printf("no generated pi/K/p found in %s\n", infile);
    return;
  }
  double plo = std::floor(pmin * 2) / 2, phi = std::ceil(pmax * 2) / 2;
  if (phi - plo < 1) {
    plo = std::max(0., plo - 0.5);
    phi += 0.5;
  }
  const double mtrue = kMass[truth];
  const double ylo   = std::max(0., std::floor(thetaExp(plo > 0 ? plo : 0.1, mtrue, n) / 20) * 20 - 40);
  const double yhi   = 1000. * std::acos(1. / n) + 20;
  const int npb      = 20;

  auto* hThP   = new TH2D("hThP", "#theta_{C} vs p;p [GeV];#theta_{C} [mrad]", 100, plo, phi, 200, ylo, yhi);
  auto* hRes   = new TH1D("hRes", TString::Format("#theta_{C} - #theta_{exp}(%s);#Delta#theta_{C} [mrad];tracks", kName[truth]), 160, -20, 20);
  auto* hResP  = new TH2D("hResP", "", npb, plo, phi, 160, -20, 20);
  auto* hNpeTh = new TH2D("hNpeTh", "photons used vs track polar angle;#theta_{track} [deg];npe", 70, 20, 160, 80, 0, 160);
  auto* pNpeTh = new TProfile("pNpeTh", "", 70, 20, 160);
  auto* hResTh = new TH2D("hResTh", "#theta_{C} - #theta_{exp} vs track polar angle;#theta_{track} [deg];#Delta#theta_{C} [mrad]", 70, 20, 160, 80, -20, 20);
  auto* pResTh = new TProfile("pResTh", "", 70, 20, 160);
  TH1D* hAll   = new TH1D("hAll", "", npb, plo, phi);
  TH1D* hId[3];
  for (int k = 0; k < 3; k++) {
    hId[k] = new TH1D(TString::Format("hId%d", k), "", npb, plo, phi);
  }

  long nGen = 0, nMatched = 0, nRec = 0, nId[3] = {0, 0, 0};
  for (Long64_t i = 0; i < nent; i++) {
    t->GetEntry(i);
    const auto st = leaf(t, "MCParticles.generatorStatus");
    const auto px = leaf(t, "MCParticles.momentum.x"), py = leaf(t, "MCParticles.momentum.y"),
               pz = leaf(t, "MCParticles.momentum.z");
    double pgen = -1;
    for (size_t k = 0; k < st.size(); k++) {
      if (st[k] == 1) {
        pgen = std::sqrt(px[k] * px[k] + py[k] * py[k] + pz[k] * pz[k]);
        break;
      }
    }
    if (pgen <= 0) {
      continue;
    }
    ++nGen;

    // segment of the generated particle: closest |p| at the bar entry
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
    ++nMatched;
    const auto k0     = static_cast<size_t>(pb[best]);
    const double thTr = std::acos(sz[k0] / pbest) * 180. / M_PI;
    const double thC  = ang[best];
    const double thE  = thetaExp(pbest, mtrue, n);
    hAll->Fill(pbest);
    if (thC <= 0) {
      continue; // no PID for this track
    }
    ++nRec;
    hThP->Fill(pbest, thC);
    hNpeTh->Fill(thTr, npe[best]);
    pNpeTh->Fill(thTr, npe[best]);
    if (thE > 0) {
      hRes->Fill(thC - thE);
      hResP->Fill(pbest, thC - thE);
      hResTh->Fill(thTr, thC - thE);
      if (std::abs(thC - thE) < 20) {
        pResTh->Fill(thTr, thC - thE);
      }
    }
    // best hypothesis (largest weight)
    const auto hb  = leaf(t, "DircGRParticleIDs.hypotheses_begin");
    const auto he  = leaf(t, "DircGRParticleIDs.hypotheses_end");
    const auto hp  = leaf(t, "_DircGRParticleIDs_hypotheses.PDG");
    const auto hw  = leaf(t, "_DircGRParticleIDs_hypotheses.weight");
    if (static_cast<size_t>(best) < hb.size()) {
      int bh   = -1;
      double w = -1e30;
      for (auto h = static_cast<size_t>(hb[best]); h < static_cast<size_t>(he[best]); h++) {
        if (hw[h] > w) {
          w  = hw[h];
          bh = species(static_cast<int>(hp[h]));
        }
      }
      if (bh >= 0) {
        hId[bh]->Fill(pbest);
        ++nId[bh];
      }
    }
  }

  // ---- page ----
  auto* c = new TCanvas("c", "DircGR QA", 1500, 950);
  c->Divide(3, 2, 0.005, 0.005);
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.045);

  // 1. thetaC vs p
  c->cd(1);
  gPad->SetRightMargin(0.12);
  gPad->SetLogz();
  hThP->Draw("colz");
  auto* leg = new TLegend(0.55, 0.15, 0.86, 0.38);
  for (int k = 0; k < 3; k++) {
    auto* fe = new TF1(TString::Format("fexp%d", k),
                       [k, n](double* x, double*) {
                         const double th = thetaExp(x[0], kMass[k], n);
                         return th > 0 ? th : -1.;
                       },
                       plo, phi, 0);
    fe->SetNpx(400);
    fe->SetLineColor(kColor[k]);
    fe->SetLineWidth(2);
    fe->SetLineStyle(k == truth ? 1 : 2);
    fe->Draw("same");
    leg->AddEntry(fe, TString::Format("%s expected (n=%.4f)", kName[k], n), "l");
  }
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  // 2. residual
  c->cd(2);
  hRes->Draw();
  double mean = 0, sig = 0;
  if (TF1* fr = fitCore(hRes)) {
    fr->SetLineColor(kRed);
    fr->Draw("same");
    mean = fr->GetParameter(1);
    sig  = std::abs(fr->GetParameter(2));
    tx.DrawLatex(0.14, 0.84, TString::Format("mean %.2f mrad", mean));
    tx.DrawLatex(0.14, 0.78, TString::Format("#sigma %.2f mrad", sig));
  }
  tx.DrawLatex(0.62, 0.84, TString::Format("%ld / %ld tracks", nRec, nMatched));

  // 3. resolution and offset vs p
  c->cd(3);
  auto* gS = new TGraphErrors();
  auto* gM = new TGraphErrors();
  for (int b = 1; b <= npb; b++) {
    TH1D* hs = hResP->ProjectionY(TString::Format("hs%d", b), b, b);
    if (TF1* fs = fitCore(hs)) {
      const double pc = hResP->GetXaxis()->GetBinCenter(b);
      gS->AddPoint(pc, std::abs(fs->GetParameter(2)));
      gS->SetPointError(gS->GetN() - 1, 0, fs->GetParError(2));
      gM->AddPoint(pc, fs->GetParameter(1));
      gM->SetPointError(gM->GetN() - 1, 0, fs->GetParError(1));
    }
  }
  auto* fr3 = gPad->DrawFrame(plo, -2, phi, 6, "track #theta_{C} resolution and offset;p [GeV];[mrad]");
  (void)fr3;
  TLine l0(plo, 0, phi, 0);
  l0.SetLineStyle(3);
  l0.DrawClone();
  gS->SetMarkerStyle(20);
  gS->SetMarkerColor(kBlack);
  gS->Draw("p");
  gM->SetMarkerStyle(24);
  gM->SetMarkerColor(kRed + 1);
  gM->SetLineColor(kRed + 1);
  gM->Draw("p");
  auto* leg3 = new TLegend(0.5, 0.75, 0.88, 0.88);
  leg3->AddEntry(gS, "#sigma (Gaussian core)", "p");
  leg3->AddEntry(gM, "mean #Delta#theta_{C}", "p");
  leg3->SetBorderSize(0);
  leg3->Draw();

  // 4. npe vs track polar angle
  c->cd(4);
  gPad->SetRightMargin(0.12);
  hNpeTh->Draw("colz");
  pNpeTh->SetMarkerStyle(20);
  pNpeTh->SetMarkerSize(0.6);
  pNpeTh->SetMarkerColor(kRed + 1);
  pNpeTh->Draw("same");

  // 5. residual vs track polar angle
  c->cd(5);
  gPad->SetRightMargin(0.12);
  hResTh->Draw("colz");
  pResTh->SetMarkerStyle(20);
  pResTh->SetMarkerSize(0.6);
  pResTh->SetMarkerColor(kRed + 1);
  pResTh->Draw("same");

  // 6. PID fractions vs p
  c->cd(6);
  gPad->DrawFrame(plo, 0, phi, 1.05, TString::Format("PID (true %s): best hypothesis;p [GeV];fraction of tracks", kName[truth]));
  auto* leg6 = new TLegend(0.6, 0.4, 0.88, 0.6);
  for (int k = 0; k < 3; k++) {
    hId[k]->Divide(hId[k], hAll, 1, 1, "B");
    hId[k]->SetMarkerStyle(20);
    hId[k]->SetMarkerColor(kColor[k]);
    hId[k]->SetLineColor(kColor[k]);
    hId[k]->Draw("e same");
    leg6->AddEntry(hId[k], TString::Format("ID as %s", kName[k]), "p");
  }
  leg6->SetBorderSize(0);
  leg6->Draw();

  c->cd(0);
  TString out = pdf;
  if (out.IsNull()) {
    out = infile;
    out.ReplaceAll(".edm4eic.root", "").ReplaceAll(".root", "");
    out += ".pdf";
  }
  c->SaveAs(out);

  printf("tracks: %ld generated, %ld matched to a bar-entry segment, %ld with thetaC\n", nGen,
         nMatched, nRec);
  printf("thetaC - expected(%s): mean %.2f mrad, sigma %.2f mrad (Gaussian core)\n", kAscii[truth],
         mean, sig);
  printf("PID: true %s; identified as pi %.1f%%, K %.1f%%, p %.1f%% of %ld matched tracks\n",
         kAscii[truth], 100. * nId[0] / std::max(1L, nMatched), 100. * nId[1] / std::max(1L, nMatched),
         100. * nId[2] / std::max(1L, nMatched), nMatched);
}

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 W.J. Llope
//
// chrcorrfit.C -- chromatic-correction fits for DircGR (EICrecon), from the
// per-photon calibration dumps written with -PDIRC:DircGR:calibFile=<file>
// (and no chrcorr loaded, so thetaC is the uncorrected angle).
//
// Same procedure as the standalone dircGR.cpp:
//   "apparent" (used by DircGR): for each zdir (1 = direct to the prism,
//      2 = via the mirror), thetaC [mrad] vs dt/L [ns/mm]; gaussian fit in
//      each dt/L slice, slices more than 2% off the expected thetaC get a
//      huge error, pol1 fit of the slice means over [INVVG_lower, INVVG_upper].
//      intercept = p0/1000 - thetaCexp/1000 [rad], slope = p1/1000.
//   "ideal" (truth diagnostics): thetaC [mrad] vs wavelength [nm]; gaussian
//      fit of the 365-375 nm slice, pol1 with the fixed slope -0.064 mrad/nm
//      through that point.
// Output: TGraphs vs track polar angle [deg], names <par>_Theta_0_0
// (pi+ 3 GeV set, DircGR chrcorrSpecies=0, chrcorrMomBin=0), plus a PDF.
//
// usage: root -l -b -q 'chrcorrfit.C("dumpdir","chrcorr_new.root","plots.pdf","old_chrcorr.root")'

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TSystem.h"
#include <cstdio>
#include <cmath>

const double cacutl      = 700;      // thetaC histogram range [mrad]
const double cacutu      = 900;
const double INVVG_lower = -0.00005; // pol1 fit range in dt/L [ns/mm]
const double INVVG_upper =  0.0001;
const double idealSlope  = -0.064;   // [mrad/nm], fixed (as in dircGR.cpp)
const double minPhotons  = 1000;     // min photon paths per zdir for a fit (dircGR.cpp: 100,
                                     // too few: forward-angle zdir -1 fits with 120-220 gave NaN/garbage)

// slices whose gaussian mean is >2% off the expected thetaC: ignore (huge error)
void CheckSlices(double expected, TH1D* h)
{
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double y = h->GetBinContent(i);
    if (y > 0 && fabs(y - expected) / expected > 0.02) h->SetBinError(i, 50.);
  }
}

// apparent fit for one zdir. Returns false if too few photons.
// Fills intercept [rad], slope, and the slice-mean histogram (for plotting).
bool FitApparent(TH2D* h, double expected, double& icpt, double& slope, TH1D*& hmean, TF1*& fpol)
{
  if (h->Integral() < minPhotons) return false;
  TF1 fG("fG", "gaus", cacutl, cacutu);
  fG.SetParameters(10, expected, 3.0);
  fG.SetParLimits(2, 0.0, 5.0);
  TObjArray* arr = new TObjArray();   // FitSlicesY makes it own the slice histograms: keep it
  h->FitSlicesY(&fG, 0, -1, 0, "Q0NRB", arr);
  hmean = (TH1D*)arr->At(1);
  hmean->SetName(Form("%s_mean", h->GetName()));
  hmean->SetDirectory(0);
  CheckSlices(expected, hmean);
  fpol = new TF1(Form("%s_pol1", h->GetName()), "pol1", INVVG_lower, INVVG_upper);
  hmean->Fit(fpol, "QNR");
  icpt  = fpol->GetParameter(0) / 1000. - expected / 1000.;
  slope = fpol->GetParameter(1) / 1000.;
  return std::isfinite(icpt) && std::isfinite(slope);
}

// ideal fit for one zdir: pol1 [mrad] vs wavelength [nm] with fixed slope
bool FitIdeal(TH2D* h, double expected, double& p0)
{
  if (h->Integral() < minPhotons) return false;
  TH1D* hp = h->ProjectionY(Form("%s_py", h->GetName()), 8, 8);   // 365-375 nm
  TF1 g("g", "gaus", 700, 900);
  g.SetParameters(10000, expected, 4);
  g.SetParLimits(2, 0, 4.5);
  hp->Fit(&g, "BQ0RN");
  double dint = g.GetParameter(1) - expected;   // dthetaC at 370 nm
  p0 = dint + idealSlope * (0. - 370.);          // pol1 intercept at 0 nm
  delete hp;
  return std::isfinite(p0);
}

void chrcorrfit(const char* dumpdir = "dump",
                const char* outname = "chrcorr.root",
                const char* pdf     = "chrcorrfit.pdf",
                const char* oldname = "")
{
  gStyle->SetOptStat(0);
  const int NANG = 15;
  const int angles[NANG] = {25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160};

  // output graphs: [zdir 0/1] apparent intercept/slope, ideal intercept/slope
  TGraph* gChrom[2][2];
  TGraph* gchr[2][2];
  for (int z = 0; z < 2; z++)
    for (int p = 0; p < 2; p++) {
      gChrom[z][p] = new TGraph();
      gChrom[z][p]->SetName(Form("gChrom%dpar%d_Theta_0_0", z + 1, p));
      gChrom[z][p]->SetTitle(Form("gChrom%dpar%d pi+ 3GeV", z + 1, p));
      gchr[z][p] = new TGraph();
      gchr[z][p]->SetName(Form("gchr%dpar%d_Theta_0_0", z + 1, p));
      gchr[z][p]->SetTitle(Form("gchr%dpar%d pi+ 3GeV", z + 1, p));
    }

  TCanvas* c = new TCanvas("c", "c", 1200, 900);
  c->Print(Form("%s[", pdf));

  printf("%6s %6s %6s | %8s %10s %9s | %8s %10s %9s | %7s %7s\n", "theta", "<rows>", "exp",
         "n1", "icpt1[mr]", "slope1", "n2", "icpt2[mr]", "slope2", "ideal1", "ideal2");

  for (int ia = 0; ia < NANG; ia++) {
    TString fn = Form("%s/S1_pi3GeV_th%d.chrcalib.root", dumpdir, angles[ia]);
    TFile* f = TFile::Open(fn);
    if (!f || f->IsZombie()) { printf("%s: missing\n", fn.Data()); continue; }
    TTree* t = (TTree*)f->Get("chrcalib");
    float thetaDeg, dtOverL, thetaC, thetaCexp, lambda;
    int zdir;
    t->SetBranchAddress("thetaDeg", &thetaDeg);
    t->SetBranchAddress("dtOverL", &dtOverL);
    t->SetBranchAddress("thetaC", &thetaC);
    t->SetBranchAddress("thetaCexp", &thetaCexp);
    t->SetBranchAddress("lambda", &lambda);
    t->SetBranchAddress("zdir", &zdir);

    TH2D* hdt[2];
    TH2D* hwl[2];
    for (int z = 0; z < 2; z++) {
      hdt[z] = new TH2D(Form("hdt%d_%d", z + 1, angles[ia]),
                        Form("#theta=%d#circ, zdir %s;(t_{hit}-t_{exp})/L [ns/mm];#theta_{C} [mrad]",
                             angles[ia], z == 0 ? "-1 (direct)" : "+1 (mirror)"),
                        44, -0.00022, 0.00022, 200, cacutl, cacutu);
      hwl[z] = new TH2D(Form("hwl%d_%d", z + 1, angles[ia]),
                        Form("#theta=%d#circ, zdir %s;#lambda [nm];#theta_{C} [mrad]",
                             angles[ia], z == 0 ? "-1 (direct)" : "+1 (mirror)"),
                        36, 295, 655, 100, cacutl, cacutu);
    }
    double sumTheta = 0, expected = 0;
    long n = t->GetEntries();
    for (long i = 0; i < n; i++) {
      t->GetEntry(i);
      int z = (zdir == -1) ? 0 : 1;
      hdt[z]->Fill(dtOverL, thetaC);
      hwl[z]->Fill(lambda, thetaC);
      sumTheta += thetaDeg;
      expected = thetaCexp;
    }
    // graph abscissa: the nominal sample angle (= the primary's polar angle,
    // what DircGR evaluates at). The row mean is printed for information only:
    // secondary tracks pull it off (e.g. 85 deg for the 80 deg sample).
    double thDeg = angles[ia];
    double thMean = (n > 0) ? sumTheta / n : 0;

    double icpt[2] = {0, 0}, slope[2] = {0, 0}, ideal[2] = {0, 0};
    bool okA[2], okI[2];
    TH1D* hmean[2] = {nullptr, nullptr};
    TF1* fpol[2]   = {nullptr, nullptr};
    for (int z = 0; z < 2; z++) {
      okA[z] = FitApparent(hdt[z], expected, icpt[z], slope[z], hmean[z], fpol[z]);
      okI[z] = FitIdeal(hwl[z], expected, ideal[z]);
      if (okA[z]) {
        gChrom[z][0]->SetPoint(gChrom[z][0]->GetN(), thDeg, icpt[z]);
        gChrom[z][1]->SetPoint(gChrom[z][1]->GetN(), thDeg, slope[z]);
      }
      if (okI[z]) {
        gchr[z][0]->SetPoint(gchr[z][0]->GetN(), thDeg, ideal[z]);
        gchr[z][1]->SetPoint(gchr[z][1]->GetN(), thDeg, idealSlope);
      }
    }
    printf("%6.1f %6.1f %6.1f | %8.0f %10.3f %9.2f | %8.0f %10.3f %9.2f | %7.2f %7.2f\n", thDeg, thMean, expected,
           hdt[0]->Integral(), 1000 * icpt[0], slope[0], hdt[1]->Integral(), 1000 * icpt[1], slope[1],
           ideal[0], ideal[1]);

    // page per angle: dt/L (top), wavelength (bottom); direct (left), mirror (right)
    c->Clear();
    c->Divide(2, 2);
    for (int z = 0; z < 2; z++) {
      c->cd(1 + z);
      hdt[z]->GetYaxis()->SetRangeUser(expected - 25, expected + 25);
      hdt[z]->Draw("colz");
      if (hmean[z]) {
        hmean[z]->SetMarkerStyle(20);
        hmean[z]->SetMarkerSize(0.6);
        hmean[z]->Draw("same e");
        fpol[z]->SetLineColor(kRed);
        fpol[z]->Draw("same");
      }
      TLine* le = new TLine(-0.00022, expected, 0.00022, expected);
      le->SetLineStyle(2);
      le->Draw();
      c->cd(3 + z);
      hwl[z]->GetYaxis()->SetRangeUser(expected - 25, expected + 25);
      hwl[z]->Draw("colz");
      if (okI[z]) {
        TF1* fi = new TF1(Form("fi%d_%d", z, angles[ia]), "pol1", 295, 655);
        fi->SetParameters(expected + ideal[z], idealSlope);
        fi->SetLineColor(kRed);
        fi->Draw("same");
      }
    }
    c->Print(pdf);
    f->Close();
  }

  // summary: new vs old, apparent intercept/slope for both zdir
  TFile* fold = (oldname && *oldname) ? TFile::Open(oldname) : nullptr; // optional comparison
  c->Clear();
  c->Divide(2, 2);
  const char* ytit[2] = {"intercept [rad]", "slope [rad / (ns/mm)]"};
  for (int z = 0; z < 2; z++)
    for (int p = 0; p < 2; p++) {
      c->cd(1 + 2 * p + z);
      TGraph* gn = gChrom[z][p];
      gn->SetMarkerStyle(20);
      gn->SetMarkerColor(kRed);
      gn->SetLineColor(kRed);
      gn->GetXaxis()->SetLimits(15, 170);
      gn->GetXaxis()->SetTitle("track #theta [deg]");
      gn->GetYaxis()->SetTitle(ytit[p]);
      gn->Draw("apl");
      TGraph* go = fold ? (TGraph*)fold->Get(gn->GetName()) : nullptr;
      double ymin = gn->GetYaxis()->GetXmin(), ymax = gn->GetYaxis()->GetXmax();
      if (go) {
        go->SetMarkerStyle(24);
        go->Draw("pl same");
        for (int i = 0; i < go->GetN(); i++) {
          ymin = std::min(ymin, go->GetY()[i]);
          ymax = std::max(ymax, go->GetY()[i]);
        }
        gn->GetYaxis()->SetRangeUser(ymin - 0.1 * (ymax - ymin), ymax + 0.1 * (ymax - ymin));
      }
      if (z == 0 && p == 0) {
        TLegend* lg = new TLegend(0.15, 0.75, 0.5, 0.88);
        lg->AddEntry(gn, "new LUTs (EICrecon dump)", "pl");
        if (go) lg->AddEntry(go, "old chrcorr.root", "pl");
        lg->Draw();
      }
    }
  c->Print(pdf);
  c->Print(Form("%s]", pdf));

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (int z = 0; z < 2; z++)
    for (int p = 0; p < 2; p++) {
      gChrom[z][p]->Write();
      gchr[z][p]->Write();
    }
  fout->Close();
  printf("wrote %s and %s\n", outname, pdf);
}

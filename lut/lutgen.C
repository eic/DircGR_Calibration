// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 W.J. Llope
//
// lutgen.C -- build one DircGR look-up-table file, LUT/lut_bar<B>.root, from
// the LUT photon-gun simulation of bar B (scripts/generateLutData.sh, files
// <indir>/bar<B>_c*.edm4hep.root).
//
// Input per MCP hit (DIRCHits, surface 0) of the gun photon: the MCP cell
// (cellID x, y), the photon direction at creation (= its direction in the
// bar, from MCParticles) and the hit time (gun point -> MCP).
//
// Path finding, per MCP pixel. From a point source the start direction fixes
// the whole trajectory, so the photons that reach one pixel along one optical
// path (one sequence of reflections in lens and prism) fill a compact patch
// in direction space, a few mrad wide; different paths are separate patches.
// The photons of a pixel are grouped by friends-of-friends (FoF): two photons
// are linked when their directions differ by less than `link` [mrad] AND
// their hit times by less than `tlink` [ns]; a group is everything connected
// through chains of links. Each group with >= `minop` photons is one path:
//   OPdir    normalised mean start direction
//   OPtime   mean hit time [ns]
//   OPprob   nOP(path) / nOP(all kept paths of the pixel)
//   OPdirRMS rms angle to the mean direction [mrad], OPtimeRMS [ns]
// Paths are stored most populated first.
//
// Pixel indices (dircLUT ipixx, ipixy) are DircGR's box-local grid, from the
// MCP cellID:  ipixx = (nPixX-1)/2 - cellY,  ipixy = cellX + (nPixY-1)/2.
// The grid TParameter<int>s nPixX, nPixY, nBar, bar are written next to the
// tree (what EICrecon's DircGRLUT requires), plus hNpath / hNop (paths and
// photons per pixel) and a TNamed with the settings.
//
// Defaults: tuned on bar 0 (10 x 2M photons): paths/pixel, fraction of photons
// kept (~97%) and path rms (~3.3 mrad median) plateau for link 1-3 mrad;
// tlink 0.05-1 ns gives the same paths (it only keeps photons that arrive ns
// apart from being merged).
//
// Compile it (ACLiC, "+"), the interpreter is much slower:
//   root -l -b -q 'lut/lutgen.C+(5,"lutsim","LUT")'
//   for b in 0 1 2 3 4 5 6 7 8 9; do root -l -b -q "lut/lutgen.C+($b,\"lutsim\",\"LUT\")"; done
// Memory ~1.1 GB and ~4 min per bar for 10 x 2M photons.

#include <TFile.h>
#include <TH2D.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace lutgen_detail {

// one field of a DD4hep cellID ("name:offset:width" or "name:width", width < 0 = signed)
struct Field {
  std::string name;
  int offset = 0, width = 0;
  bool isSigned = false;
  long long get(ULong64_t id) const {
    const ULong64_t mask = (width >= 64) ? ~0ULL : ((1ULL << width) - 1);
    long long v          = (long long)((id >> offset) & mask);
    if (isSigned && v >= (1LL << (width - 1))) v -= (1LL << width);
    return v;
  }
};

std::vector<Field> parseEncoding(const std::string& enc) {
  std::vector<Field> out;
  int next = 0;
  size_t p = 0;
  while (p < enc.size()) {
    size_t q = enc.find(',', p);
    if (q == std::string::npos) q = enc.size();
    const std::string tok = enc.substr(p, q - p);
    p                     = q + 1;
    std::vector<std::string> parts;
    size_t a = 0;
    while (true) {
      size_t b = tok.find(':', a);
      parts.push_back(tok.substr(a, b == std::string::npos ? std::string::npos : b - a));
      if (b == std::string::npos) break;
      a = b + 1;
    }
    Field f;
    f.name = parts[0];
    int w  = 0;
    if (parts.size() == 3) {
      f.offset = std::stoi(parts[1]);
      w        = std::stoi(parts[2]);
    } else if (parts.size() == 2) {
      f.offset = next;
      w        = std::stoi(parts[1]);
    } else {
      continue;
    }
    f.isSigned = w < 0;
    f.width    = std::abs(w);
    next       = f.offset + f.width;
    out.push_back(f);
  }
  return out;
}

// DIRCHits cellID encoding from the podio "metadata" tree (written by npsim)
std::string readEncoding(const char* fname) {
  TFile f(fname);
  auto* t = f.Get<TTree>("metadata");
  if (t == nullptr) return "";
  std::vector<std::string>* keys               = nullptr;
  std::vector<std::vector<std::string>>* value = nullptr;
  if (t->GetBranch("GPStringKeys") == nullptr || t->GetBranch("GPStringValues") == nullptr)
    return "";
  t->SetBranchAddress("GPStringKeys", &keys);
  t->SetBranchAddress("GPStringValues", &value);
  std::string enc;
  if (t->GetEntry(0) > 0) {
    for (size_t i = 0; i < keys->size(); i++) {
      if ((*keys)[i].rfind("DIRCHits", 0) == 0 &&
          (*keys)[i].find("CellIDEncoding") != std::string::npos && !(*value)[i].empty())
        enc = (*value)[i][0];
    }
  }
  t->ResetBranchAddresses();
  return enc;
}

struct Hit {
  int pix;
  float dx, dy, dz, t;
};

// union-find
int findRoot(std::vector<int>& par, int i) {
  while (par[i] != i) {
    par[i] = par[par[i]];
    i      = par[i];
  }
  return i;
}
void unite(std::vector<int>& par, int a, int b) {
  a = findRoot(par, a);
  b = findRoot(par, b);
  if (a != b) par[std::max(a, b)] = std::min(a, b);
}

} // namespace lutgen_detail

int lutgen(int bar = 0, const char* indir = "lutsim", const char* outdir = "LUT",
           double link = 2.0, double tlink = 0.2, int minop = 20, int nPixX = 107,
           int nPixY = 73, int nBar = 10, int maxfiles = 0) {
  using namespace lutgen_detail;

  // ---- input files, sorted
  std::vector<std::string> files;
  {
    TSystemDirectory dir(indir, indir);
    const TString pre = TString::Format("bar%d_c", bar);
    if (auto* l = dir.GetListOfFiles()) {
      for (auto* o : *l) {
        const TString n = o->GetName();
        if (n.BeginsWith(pre) && n.EndsWith(".edm4hep.root"))
          files.push_back(std::string(indir) + "/" + n.Data());
      }
      delete l;
    }
  }
  std::sort(files.begin(), files.end());
  if (maxfiles > 0 && (int)files.size() > maxfiles) files.resize(maxfiles);
  if (files.empty()) {
    printf("lutgen: no %s/bar%d_c*.edm4hep.root\n", indir, bar);
    return 1;
  }
  if (nPixX % 2 == 0 || nPixY % 2 == 0) {
    printf("lutgen: nPixX, nPixY must be odd (MCP grid centred on cell 0)\n");
    return 1;
  }
  const int halfX = (nPixX - 1) / 2; // cell y range -halfX..halfX
  const int halfY = (nPixY - 1) / 2; // cell x range -halfY..halfY

  // ---- cellID fields
  std::string enc = readEncoding(files[0].c_str());
  if (enc.empty()) {
    enc = "system:8,module:4,surface:1,section:2,bar:4,x:-16,y:-16";
    printf("lutgen: no DIRCHits CellIDEncoding in %s, using %s\n", files[0].c_str(), enc.c_str());
  }
  const auto fields = parseEncoding(enc);
  auto fieldOf      = [&](const char* nm) -> const Field* {
    for (auto& f : fields)
      if (f.name == nm) return &f;
    return nullptr;
  };
  const Field *fMod = fieldOf("module"), *fSurf = fieldOf("surface"), *fX = fieldOf("x"),
              *fY = fieldOf("y");
  if (!fMod || !fSurf || !fX || !fY) {
    printf("lutgen: cellID encoding '%s' lacks module/surface/x/y\n", enc.c_str());
    return 1;
  }

  // ---- read the hits
  std::vector<Hit> hits;
  hits.reserve(16000000);
  long long nthrown = 0, nbadprim = 0, nout = 0;
  double vx0 = 0, vy0 = 0, vz0 = 0, maxdv = 0;
  for (const auto& fn : files) {
    TFile f(fn.c_str());
    auto* t = f.Get<TTree>("events");
    if (t == nullptr) {
      printf("lutgen: no events tree in %s\n", fn.c_str());
      return 1;
    }
    TTreeReader r(t);
    TTreeReaderArray<ULong_t> cid(r, "DIRCHits.cellID");
    TTreeReaderArray<Float_t> htime(r, "DIRCHits.time");
    TTreeReaderArray<Int_t> hlink(r, "_DIRCHits_particle.index");
    TTreeReaderArray<Int_t> gst(r, "MCParticles.generatorStatus");
    TTreeReaderArray<Double_t> px(r, "MCParticles.momentum.x"), py(r, "MCParticles.momentum.y"),
        pz(r, "MCParticles.momentum.z");
    TTreeReaderArray<Double_t> vx(r, "MCParticles.vertex.x"), vy(r, "MCParticles.vertex.y"),
        vz(r, "MCParticles.vertex.z");
    while (r.Next()) {
      nthrown++;
      int ip = -1, nprim = 0;
      for (size_t i = 0; i < gst.GetSize(); i++)
        if (gst[i] == 1) {
          ip = (int)i;
          nprim++;
        }
      if (nprim != 1) {
        nbadprim++;
        continue;
      }
      if (nthrown == 1) {
        vx0 = vx[ip];
        vy0 = vy[ip];
        vz0 = vz[ip];
      }
      maxdv = std::max(maxdv, std::hypot(vx[ip] - vx0, vy[ip] - vy0, vz[ip] - vz0));
      const double p = std::sqrt(px[ip] * px[ip] + py[ip] * py[ip] + pz[ip] * pz[ip]);
      const float dx = px[ip] / p, dy = py[ip] / p, dz = pz[ip] / p;
      for (size_t k = 0; k < cid.GetSize(); k++) {
        if (hlink[k] != ip) continue;
        const ULong64_t c = cid[k];
        if (fSurf->get(c) != 0) continue; // MCP surface only
        if (fMod->get(c) != 0) {
          printf("lutgen: MCP hit in bar box %lld: the LUT gun must be in box 0\n", fMod->get(c));
          return 1;
        }
        const int ix = halfX - (int)fY->get(c);
        const int iy = (int)fX->get(c) + halfY;
        if (ix < 0 || ix >= nPixX || iy < 0 || iy >= nPixY) {
          nout++;
          continue;
        }
        hits.push_back({ix * nPixY + iy, dx, dy, dz, htime[k]});
      }
    }
  }
  if (nout > 0) {
    printf("lutgen: %lld hits outside the %d x %d grid: wrong nPixX/nPixY?\n", nout, nPixX, nPixY);
    return 1;
  }
  if (maxdv > 1e-6)
    printf("lutgen: WARNING gun position varies by %.3g mm between events\n", maxdv);
  printf("bar %d: %zu files, thrown %lld, MCP hits %zu (%.1f%%), gun at (%.3f, %.3f, %.3f) mm, "
         "%lld events without a single primary\n",
         bar, files.size(), nthrown, hits.size(), 100. * hits.size() / std::max(1LL, nthrown), vx0,
         vy0, vz0, nbadprim);

  // ---- group by pixel, keeping the read order (counting sort)
  const int npix = nPixX * nPixY;
  std::vector<int> pbeg(npix + 1, 0);
  for (auto& h : hits)
    pbeg[h.pix + 1]++;
  for (int i = 0; i < npix; i++)
    pbeg[i + 1] += pbeg[i];
  std::vector<int> order(hits.size());
  {
    std::vector<int> fill(pbeg.begin(), pbeg.end() - 1);
    for (size_t i = 0; i < hits.size(); i++)
      order[fill[hits[i].pix]++] = (int)i;
  }

  // ---- output
  gSystem->mkdir(outdir, true);
  const TString fout = TString::Format("%s/lut_bar%d.root", outdir, bar);
  TFile out(fout, "RECREATE");
  if (out.IsZombie()) {
    printf("lutgen: cannot write %s\n", fout.Data());
    return 1;
  }
  TTree tree("dircLUT", "DircGR LUT: optical-photon paths per MCP pixel");
  int T_ipixx, T_ipixy, T_npath, T_ipath, T_nOP;
  double T_prob, T_dx, T_dy, T_dz, T_t, T_drms, T_trms;
  tree.Branch("ipixx", &T_ipixx, "ipixx/I");
  tree.Branch("ipixy", &T_ipixy, "ipixy/I");
  tree.Branch("npath", &T_npath, "npath/I");
  tree.Branch("ipath", &T_ipath, "ipath/I");
  tree.Branch("nOP", &T_nOP, "nOP/I");
  tree.Branch("OPprob", &T_prob, "OPprob/D");
  tree.Branch("OPdirx", &T_dx, "OPdirx/D");
  tree.Branch("OPdiry", &T_dy, "OPdiry/D");
  tree.Branch("OPdirz", &T_dz, "OPdirz/D");
  tree.Branch("OPtime", &T_t, "OPtime/D");
  tree.Branch("OPdirRMS", &T_drms, "OPdirRMS/D");
  tree.Branch("OPtimeRMS", &T_trms, "OPtimeRMS/D");
  TH2D hNpath("hNpath", "paths per pixel;ipixx;ipixy", nPixX, -0.5, nPixX - 0.5, nPixY, -0.5,
              nPixY - 0.5);
  TH2D hNop("hNop", "photons per pixel;ipixx;ipixy", nPixX, -0.5, nPixX - 0.5, nPixY, -0.5,
            nPixY - 0.5);

  // ---- FoF per pixel
  const double L   = link * 1e-3; // chord ~ angle [rad] for small angles
  const double L2  = L * L;
  long long nkeptTot = 0, npathTot = 0, npixHit = 0;
  int npathMax       = 0;
  std::vector<double> rmsAll;
  for (int p = 0; p < npix; p++) {
    const int n = pbeg[p + 1] - pbeg[p];
    if (n == 0) continue;
    npixHit++;
    const int* idx = &order[pbeg[p]];
    std::vector<double> d(3 * n), t(n);
    for (int i = 0; i < n; i++) {
      const Hit& h  = hits[idx[i]];
      const double x = h.dx, y = h.dy, z = h.dz, nn = std::sqrt(x * x + y * y + z * z);
      d[3 * i]     = x / nn;
      d[3 * i + 1] = y / nn;
      d[3 * i + 2] = z / nn;
      t[i]         = h.t;
    }
    // cells of size L in direction space: linked photons are in the same or
    // neighbouring cells, so only those pairs are tested
    std::vector<std::pair<int64_t, int>> cell(n);
    auto ci = [&](double v) { return (int64_t)std::floor(v / L) + (1 << 20); };
    auto key = [](int64_t a, int64_t b, int64_t c) { return (a << 42) | (b << 21) | c; };
    for (int i = 0; i < n; i++)
      cell[i] = {key(ci(d[3 * i]), ci(d[3 * i + 1]), ci(d[3 * i + 2])), i};
    std::sort(cell.begin(), cell.end());
    std::vector<int> par(n);
    std::iota(par.begin(), par.end(), 0);
    auto linked = [&](int i, int j) {
      const double ax = d[3 * i] - d[3 * j], ay = d[3 * i + 1] - d[3 * j + 1],
                   az = d[3 * i + 2] - d[3 * j + 2];
      return ax * ax + ay * ay + az * az <= L2 && std::fabs(t[i] - t[j]) < tlink;
    };
    for (size_t a = 0; a < cell.size();) {
      size_t b = a;
      while (b < cell.size() && cell[b].first == cell[a].first)
        b++;
      const int64_t k  = cell[a].first;
      const int64_t cx = (k >> 42) & 0x1FFFFF, cy = (k >> 21) & 0x1FFFFF, cz = k & 0x1FFFFF;
      // pairs inside the cell
      for (size_t i = a; i < b; i++)
        for (size_t j = i + 1; j < b; j++)
          if (linked(cell[i].second, cell[j].second)) unite(par, cell[i].second, cell[j].second);
      // pairs with the 13 "forward" neighbour cells (each cell pair once)
      for (int ox = -1; ox <= 1; ox++)
        for (int oy = -1; oy <= 1; oy++)
          for (int oz = -1; oz <= 1; oz++) {
            const int64_t nk = key(cx + ox, cy + oy, cz + oz);
            if (nk <= k) continue;
            auto lo = std::lower_bound(cell.begin(), cell.end(), std::make_pair(nk, -1));
            for (auto it = lo; it != cell.end() && it->first == nk; ++it)
              for (size_t i = a; i < b; i++)
                if (linked(cell[i].second, it->second)) unite(par, cell[i].second, it->second);
          }
      a = b;
    }
    // groups, labelled in order of their first photon
    std::vector<int> lab(n, -1), rootLab(n, -1), size;
    for (int i = 0; i < n; i++) {
      const int r = findRoot(par, i);
      if (rootLab[r] < 0) {
        rootLab[r] = (int)size.size();
        size.push_back(0);
      }
      lab[i] = rootLab[r];
      size[lab[i]]++;
    }
    std::vector<int> keep;
    int nkept = 0;
    for (int g = 0; g < (int)size.size(); g++)
      if (size[g] >= minop) {
        keep.push_back(g);
        nkept += size[g];
      }
    std::stable_sort(keep.begin(), keep.end(), [&](int a, int b) { return size[a] > size[b]; });
    T_ipixx = p / nPixY;
    T_ipixy = p % nPixY;
    T_npath = (int)keep.size();
    hNpath.SetBinContent(T_ipixx + 1, T_ipixy + 1, T_npath);
    hNop.SetBinContent(T_ipixx + 1, T_ipixy + 1, n);
    nkeptTot += nkept;
    npathTot += T_npath;
    npathMax = std::max(npathMax, T_npath);
    for (int ip = 0; ip < T_npath; ip++) {
      const int g = keep[ip];
      double mx = 0, my = 0, mz = 0, mt = 0;
      for (int i = 0; i < n; i++)
        if (lab[i] == g) {
          mx += d[3 * i];
          my += d[3 * i + 1];
          mz += d[3 * i + 2];
          mt += t[i];
        }
      const double mn = std::sqrt(mx * mx + my * my + mz * mz);
      mx /= mn;
      my /= mn;
      mz /= mn;
      mt /= size[g];
      double sa = 0, st = 0;
      for (int i = 0; i < n; i++)
        if (lab[i] == g) {
          const double c =
              std::clamp(d[3 * i] * mx + d[3 * i + 1] * my + d[3 * i + 2] * mz, -1.0, 1.0);
          const double ang = std::acos(c);
          sa += ang * ang;
          st += (t[i] - mt) * (t[i] - mt);
        }
      T_ipath = ip;
      T_nOP   = size[g];
      T_prob  = (double)size[g] / nkept;
      T_dx    = mx;
      T_dy    = my;
      T_dz    = mz;
      T_t     = mt;
      T_drms  = 1e3 * std::sqrt(sa / size[g]);
      T_trms  = std::sqrt(st / size[g]);
      rmsAll.push_back(T_drms);
      tree.Fill();
    }
  }

  out.cd();
  tree.Write();
  hNpath.Write();
  hNop.Write();
  TParameter<int>("nPixX", nPixX).Write("nPixX", TObject::kOverwrite);
  TParameter<int>("nPixY", nPixY).Write("nPixY", TObject::kOverwrite);
  TParameter<int>("nBar", nBar).Write("nBar", TObject::kOverwrite);
  TParameter<int>("bar", bar).Write("bar", TObject::kOverwrite);
  TNamed("lutgen_info",
         TString::Format("lutgen.C: bar %d, %zu files from %s, thrown %lld, MCP hits %zu, gun "
                         "(%.3f,%.3f,%.3f) mm, link %g mrad, tlink %g ns, minop %d",
                         bar, files.size(), indir, nthrown, hits.size(), vx0, vy0, vz0, link,
                         tlink, minop)
             .Data())
      .Write();
  std::nth_element(rmsAll.begin(), rmsAll.begin() + rmsAll.size() / 2, rmsAll.end());
  printf("wrote %s: %lld pixels, %lld paths (paths/pixel mean %.1f, max %d), photons in paths "
         "%.1f%%, path rms angle median %.2f mrad\n",
         fout.Data(), npixHit, npathTot, (double)npathTot / std::max(1LL, npixHit), npathMax,
         100. * nkeptTot / std::max<size_t>(1, hits.size()),
         rmsAll.empty() ? 0. : rmsAll[rmsAll.size() / 2]);
  out.Close();
  return 0;
}

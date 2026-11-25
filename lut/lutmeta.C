// lutmeta.C: write the DircGR LUT grid parameters into lut_bar{0..nBar-1}.root
//
// Each file gets TParameter<int> nPixX, nPixY, nBar (the same in all files)
// and bar (= the file index). EICrecon's DircGRLUT reads the grid from these
// (the dircLUT tree itself is not touched). Existing parameters are
// overwritten.
//
// Grid = DIRCHits MCP cellID grid: 107 (cell y) x 73 (cell x), 10 bars.
//
//   root -l -b -q 'lutmeta.C("LUT")'
//   root -l -b -q 'lutmeta.C("LUT", 107, 73, 10)'

#include <TFile.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>

#include <cstdio>

int lutmeta(const char* dir = "LUT", int nPixX = 107, int nPixY = 73, int nBar = 10) {
  int nbad = 0;
  for (int bar = 0; bar < nBar; bar++) {
    const TString fname = TString::Format("%s/lut_bar%d.root", dir, bar);
    TFile f(fname, "UPDATE");
    if (f.IsZombie()) {
      printf("lutmeta: cannot open %s\n", fname.Data());
      nbad++;
      continue;
    }
    auto* t = f.Get<TTree>("dircLUT");
    if (t == nullptr) {
      printf("lutmeta: no dircLUT tree in %s\n", fname.Data());
      nbad++;
      continue;
    }
    // check the pixels fit the grid
    const double xmax = t->GetMaximum("ipixx"), ymax = t->GetMaximum("ipixy");
    const double xmin = t->GetMinimum("ipixx"), ymin = t->GetMinimum("ipixy");
    if (xmin < 0 || ymin < 0 || xmax >= nPixX || ymax >= nPixY) {
      printf("lutmeta: %s pixels x %g..%g y %g..%g outside %dx%d\n", fname.Data(), xmin, xmax,
             ymin, ymax, nPixX, nPixY);
      nbad++;
      continue;
    }
    TParameter<int>("nPixX", nPixX).Write("nPixX", TObject::kOverwrite);
    TParameter<int>("nPixY", nPixY).Write("nPixY", TObject::kOverwrite);
    TParameter<int>("nBar", nBar).Write("nBar", TObject::kOverwrite);
    TParameter<int>("bar", bar).Write("bar", TObject::kOverwrite);
    printf("%s: %lld paths, pixels x %g..%g y %g..%g; wrote nPixX=%d nPixY=%d nBar=%d bar=%d\n",
           fname.Data(), t->GetEntries(), xmin, xmax, ymin, ymax, nPixX, nPixY, nBar, bar);
  }
  return nbad;
}

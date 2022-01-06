// Make a simple spectrum plot
// cafe plotVariable.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"

#include "StandardRecord/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"

using namespace ana;

void plotVariable()
{
  std::cout << "\033[31m" << "BEGIN " << "\033[0m"  << std::endl;

  // File from CAkMaker dunetpc v07_09_00
  const std::string fnameNonSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nu.root";
  const std::string fnameNueSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nue.root";
  const std::string fnameTauSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nutau.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);

  const Var kRecoEnergy = SIMPLEVAR(Ev_reco_numu);
  // Variable written into CAF file for testing
  const Var kSelTrackPandizzleScore = SIMPLEVAR(selTrackPandizzleScore);


  std::cout << "\033[31m" << "AFTER JAM VARIABLE DEFINED" << "\033[0m"  << std::endl;

  // Creating variable distribution - apparently have to apply cut so pick one that passes test event
  const Binning binsJam = Binning::Simple(50, -2, 2);
  const HistAxis axJam("Jam", binsJam, kSelTrackPandizzleScore);
  Spectrum sJam(loader, axJam, !kIsAntiNu);

  std::cout << "\033[31m" << "AFTER JAM SPCTRUM CREATED" << "\033[0m"  << std::endl;

  // This is the call that actually fills in those spectra
  loader.Go();

  std::cout << "\033[31m" << "AFTER SPECTRUM LOADED" << "\033[0m"  << std::endl;

  // POT/yr * 3.5yrs * mass correction for the workspace geometry
  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // For plotting purposes we can convert to TH1s
  sJam.ToTH1(pot)->Draw("hist");

}

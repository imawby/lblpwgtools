// Make oscillated predictions
// cafe demo1.C

#include "CAFAna/Core/SpectrumLoader.h"
//#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
//#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "StandardRecord/SRProxy.h"
#include "TCanvas.h"
#include "TH1.h"

// New includes for this macro
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/OscCalcPMNSOpt.h"

using namespace ana;

void demo1()
{
  // See demo0.C for explanation of these repeated parts
  //const std::string fnameNonSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nu.root";
  //const std::string fnameNueSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nue.root";
  //const std::string fnameTauSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nutau.root";
  //SpectrumLoader loaderNonSwap(fnameNonSwap);
  //SpectrumLoader loaderNueSwap(fnameNueSwap);
  //SpectrumLoader loaderTauSwap(fnameTauSwap);
  const std::string full = "/dune/app/users/imawby/selection/fcl/MCC11/caf_hist.root";
  SpectrumLoader loaderFull(full);

  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);
  const Var kIntPdg = SIMPLEVAR(nuPDG);
  const Var kIsCC = SIMPLEVAR(isCC);
  const Var kTrackPandizzle = SIMPLEVAR(selTrackPandizzleScore);
  const Var kShowerPandrizzle = SIMPLEVAR(selShowerPandrizzleScore);

  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axEnergy("Reco energy (GeV)", binsEnergy, kRecoNumuEnergy);
  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // A cut is structured like a Var, but returning bool
  const Cut kPassesCVN([](const caf::SRProxy* sr)
                       {
                         return sr->cvnnumu > 0.5;
                       });

  // In many cases it's easier to form them from existing Vars like this
  //  const Cut kPassesCVN = kCVNNumu > 0.5;

  // A Prediction is an object holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  /*
  PredictionNoExtrap pred(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axEnergy, kIsNumuSelected);
  */

  Spectrum pred(loaderFull, axEnergy, kIsNumuSelected);

  // These calls will fill all of the constituent parts of the prediction
  /*
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();
  */
  loaderFull.Go();

  // We can extract a total prediction unoscillated
  //const Spectrum sUnosc = pred.PredictUnoscillated();
  // Or oscillated, in this case using reasonable parameters from
  // Analysis/Calcs.h
  osc::IOscCalcAdjustable* calc = DefaultOscCalc();
  //calc->SetdCP(0);
  //const Spectrum sOsc = pred.Predict(calc);
  
/*
  const Spectrum sOsc = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kBoth,
                                                  Sign::kBoth);
  

  // And we can break things down by flavour
  
  const Spectrum sUnoscNC = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);
  
  // Plot what we have so far
  //sUnosc.ToTH1(pot)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.ToTH1(pot, kRed)->Draw("hist same");
*/

pred.ToTH1(pot, kBlue)->Draw("hist same");
  // "Fake" data is synonymous with the Asimov data sample
  /*
  new TCanvas;
  sOsc.ToTH1(pot, kRed)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.FakeData(pot).ToTH1(pot)->Draw("ep same");

  // While "mock" data has statistical fluctuations in
  new TCanvas;
  sOsc.ToTH1(pot, kRed)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.MockData(pot).ToTH1(pot)->Draw("ep same");
  */
}

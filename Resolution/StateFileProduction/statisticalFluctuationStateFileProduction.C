// Make resolution plot <- with just statistical fluctuations
// cafe deltaCPResolution_statisticalOnly.C

#include "CAFAna/Analysis/Resolution.h"

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "StandardRecord/SRProxy.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMath.h"

// New includes for this macro
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/MultiExperiment.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

// also add in some RHC files!

using namespace ana;

const std::string outputFileName = "Test.root";

void statisticalFluctuationStateFileProduction();

void statisticalFluctuationStateFileProduction()
{
    std::cout << "ISOBEL YOU NEED TO IMPLEMENT THE FIDUCIAL CUT TOO" << std::endl;

   std::cout << "Reading caf files..." << std::endl;

  const std::string fnameNonSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nu.root";
  const std::string fnameNueSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nue.root";
  const std::string fnameTauSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nutau.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);

  std::cout << "Read files!" << std::endl;

  const Var kRecoNueEnergy = SIMPLEVAR(Ev_reco_nue);
  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);

  const Var kTrackPandizzle = SIMPLEVAR(selTrackPandizzleScore);
  const Var kShowerPandrizzle = SIMPLEVAR(selShowerPandrizzleScore);

  const Cut kPassesCVN_nue = SIMPLEVAR(cvnnue) > 0.85;
  const Cut kPassesCVN_numu = SIMPLEVAR(cvnnumu) > 0.5;

  // TODO - change the binning
  const Binning binsEnergy = Binning::Simple(80, 0, 10);
  const HistAxis axNueEnergy("Reco #nu_{e} energy (GeV)", binsEnergy, kRecoNueEnergy);
  const HistAxis axNumuEnergy("Reco #nu_{#mu} energy (GeV)", binsEnergy, kRecoNumuEnergy);

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap predNue_FHC_IZZLE(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kIsNueSelected);
  PredictionNoExtrap predNumu_FHC_IZZLE(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kIsNumuSelected);
  PredictionNoExtrap predNue_FHC_CVN(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kPassesCVN_nue);
  PredictionNoExtrap predNumu_FHC_CVN(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNumuEnergy, kPassesCVN_numu);

  std::cout << "Filling spectra..." << std::endl;

  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();

  std::cout << "Spectra filled!" << std::endl;

  std::cout << "Saving to file: " << outputFileName << std::endl;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  predNue_FHC_IZZLE.SaveTo(outputFile, "predNue_FHC_IZZLE");
  predNumu_FHC_IZZLE.SaveTo(outputFile, "predNumu_FHC_IZZLE");
  predNue_FHC_CVN.SaveTo(outputFile, "predNue_FHC_CVN");
  predNumu_FHC_CVN.SaveTo(outputFile, "predNumu_FHC_CVN");

  std::cout << "All done making state files..." << std::endl;
}

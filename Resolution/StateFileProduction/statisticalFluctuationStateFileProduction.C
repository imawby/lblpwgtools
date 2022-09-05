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

const std::string outputFileName = "/dune/data/users/imawby/standardCAF/StateFiles.root";

void statisticalFluctuationStateFileProduction();

void statisticalFluctuationStateFileProduction()
{
   std::cout << "Reading caf files..." << std::endl;
   
   const std::string fnameNonSwapFHC = "/dune/data/users/imawby/standardCAF/nu/caf_hist_nu.root";
   const std::string fnameNueSwapFHC = "/dune/data/users/imawby/standardCAF/nue/caf_hist_nue.root";
   const std::string fnameTauSwapFHC = "/dune/data/users/imawby/standardCAF/nutau/caf_hist_nutau.root";
   const std::string fnameNonSwapRHC = "/dune/data/users/imawby/standardCAF/anu/caf_hist_anu.root";
   const std::string fnameNueSwapRHC = "/dune/data/users/imawby/standardCAF/anue/caf_hist_anue.root";
   const std::string fnameTauSwapRHC = "/dune/data/users/imawby/standardCAF/anutau/caf_hist_anutau.root";

   /*
   const std::string fnameNonSwapFHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_FHC_nonswap.root";
   const std::string fnameNueSwapFHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_FHC_nueswap.root";
   const std::string fnameTauSwapFHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_FHC_tauswap.root";

   const std::string fnameNonSwapRHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_RHC_nonswap.root";
   const std::string fnameNueSwapRHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_RHC_nueswap.root";
   const std::string fnameTauSwapRHC = "/dune/data/users/marshalc/CAFs/obsolete/mcc11_v3/FD_RHC_tauswap.root";
   */

  SpectrumLoader loaderNonSwapFHC(fnameNonSwapFHC);
  SpectrumLoader loaderNueSwapFHC(fnameNueSwapFHC);
  SpectrumLoader loaderTauSwapFHC(fnameTauSwapFHC);

  SpectrumLoader loaderNonSwapRHC(fnameNonSwapRHC);
  SpectrumLoader loaderNueSwapRHC(fnameNueSwapRHC);
  SpectrumLoader loaderTauSwapRHC(fnameTauSwapRHC);

  std::cout << "Read files!" << std::endl;

  const Var kRecoNueEnergy = SIMPLEVAR(Ev_reco_nue);
  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);

  // TODO - change the binning
  const Binning binsEnergy = Binning::Simple(32, 0, 8);
  const HistAxis axNueEnergy("Reco #nu_{e} energy (GeV)", binsEnergy, kRecoNueEnergy);
  const HistAxis axNumuEnergy("Reco #nu_{#mu} energy (GeV)", binsEnergy, kRecoNumuEnergy);

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.

  PredictionNoExtrap predNue_FHC_IZZLE(loaderNonSwapFHC, loaderNueSwapFHC, loaderTauSwapFHC, axNueEnergy, kIsNueSelectedFHC);
  PredictionNoExtrap predNumu_FHC_IZZLE(loaderNonSwapFHC, loaderNueSwapFHC, loaderTauSwapFHC, axNumuEnergy, kIsNumuSelectedFHC);
  PredictionNoExtrap predNue_RHC_IZZLE(loaderNonSwapRHC, loaderNueSwapRHC, loaderTauSwapRHC, axNueEnergy, kIsNueSelectedRHC);
  PredictionNoExtrap predNumu_RHC_IZZLE(loaderNonSwapRHC, loaderNueSwapRHC, loaderTauSwapRHC, axNumuEnergy, kIsNumuSelectedRHC);

  PredictionNoExtrap predNue_FHC_CVN(loaderNonSwapFHC, loaderNueSwapFHC, loaderTauSwapFHC, axNueEnergy, kPassFD_CVN_NUE);
  PredictionNoExtrap predNumu_FHC_CVN(loaderNonSwapFHC, loaderNueSwapFHC, loaderTauSwapFHC, axNumuEnergy, kPassFD_CVN_NUMU);
  PredictionNoExtrap predNue_RHC_CVN(loaderNonSwapRHC, loaderNueSwapRHC, loaderTauSwapRHC, axNueEnergy, kPassFD_CVN_NUE);
  PredictionNoExtrap predNumu_RHC_CVN(loaderNonSwapRHC, loaderNueSwapRHC, loaderTauSwapRHC, axNumuEnergy, kPassFD_CVN_NUMU);

  std::cout << "Filling spectra..." << std::endl;

  loaderNonSwapFHC.Go();
  loaderNueSwapFHC.Go();
  loaderTauSwapFHC.Go();

  loaderNonSwapRHC.Go();
  loaderNueSwapRHC.Go();
  loaderTauSwapRHC.Go();

  std::cout << "Spectra filled!" << std::endl;

  std::cout << "Saving to file: " << outputFileName << std::endl;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");

  predNue_FHC_IZZLE.SaveTo(outputFile, "predNue_FHC_IZZLE");
  predNumu_FHC_IZZLE.SaveTo(outputFile, "predNumu_FHC_IZZLE");
  predNue_RHC_IZZLE.SaveTo(outputFile, "predNue_RHC_IZZLE");
  predNumu_RHC_IZZLE.SaveTo(outputFile, "predNumu_RHC_IZZLE");

  predNue_FHC_CVN.SaveTo(outputFile, "predNue_FHC_CVN");
  predNumu_FHC_CVN.SaveTo(outputFile, "predNumu_FHC_CVN");
  predNue_RHC_CVN.SaveTo(outputFile, "predNue_RHC_CVN");
  predNumu_RHC_CVN.SaveTo(outputFile, "predNumu_RHC_CVN");

  std::cout << "All done making state files..." << std::endl;
}

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
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Systs/AnaSysts.h"
#include "CAFAna/Core/Loaders.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

// also add in some RHC files!

using namespace ana;

const std::string outputFileName = "/dune/data/users/imawby/standardCAF/StateFilesDetectorSystematics.root";

void statisticalFluctuationStateFileProductionSystematics();

void statisticalFluctuationStateFileProductionSystematics()
{
   std::cout << "Reading caf files..." << std::endl;
   
  const std::string fnameNonSwapFHC = "/dune/data/users/imawby/standardCAF/nu/caf_hist_nu.root";
  const std::string fnameNueSwapFHC = "/dune/data/users/imawby/standardCAF/nue/caf_hist_nue.root";
  const std::string fnameTauSwapFHC = "/dune/data/users/imawby/standardCAF/nutau/caf_hist_nutau.root";

  const std::string fnameNonSwapRHC = "/dune/data/users/imawby/standardCAF/anu/caf_hist_anu.root";
  const std::string fnameNueSwapRHC = "/dune/data/users/imawby/standardCAF/anue/caf_hist_anue.root";
  const std::string fnameTauSwapRHC = "/dune/data/users/imawby/standardCAF/anutau/caf_hist_anutau.root";

  Loaders loadersFHC;
  loadersFHC.SetLoaderPath(fnameNonSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNonSwap);
  loadersFHC.SetLoaderPath(fnameNueSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNueSwap);
  loadersFHC.SetLoaderPath(fnameTauSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNuTauSwap);

  /*
  Loaders loadersRHC;
  loadersRHC.SetLoaderPath(fnameNonSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNonSwap);
  loadersRHC.SetLoaderPath(fnameNueSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNueSwap);
  loadersRHC.SetLoaderPath(fnameTauSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNuTauSwap);
  */
  const Var kRecoNueEnergy = SIMPLEVAR(Ev_reco_nue);
  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);

  // TODO - change the binning
  const Binning binsEnergy = Binning::Simple(32, 0, 8);
  const HistAxis axNueEnergy("Reco #nu_{e} energy (GeV)", binsEnergy, kRecoNueEnergy);
  const HistAxis axNumuEnergy("Reco #nu_{#mu} energy (GeV)", binsEnergy, kRecoNumuEnergy);

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  /*
  std::vector<const ISyst *> GetListOfSysts(bool fluxsyst_Nov17, bool xsecsyst,
                                            bool detsyst, bool useND, bool useFD,
                                            bool useNueOnE, bool useFakeDataDials,
                                            bool fluxsyst_CDR, int NFluxSysts,
                                            bool removeFDNonFitDials) 
  */

  // Only detector systematics
  std::vector<ana::ISyst const *> systematicsVector = GetListOfSysts(false, false, true, false, true, false, false, false);

  for (ana::ISyst const * sys : systematicsVector)
  {
      std::cout << "ShortName: " << sys->ShortName() << std::endl;
      std::cout << "CentralValue: " << sys->Central() << std::endl;
      std::cout << "Min: " << sys->Min() << std::endl;
      std::cout << "Max: " << sys->Max() << std::endl;
      std::cout << "Apply shift: " << sys->ApplyPenalty() << std::endl;
  }

  osc::IOscCalcAdjustable* calc = DefaultOscCalc();
  NoExtrapPredictionGenerator genNue_FHC_IZZLE(axNueEnergy, kIsNueSelectedFHC, kUnweighted);
  //NoExtrapPredictionGenerator genNue_RHC_IZZLE(axNueEnergy, kIsNueSelectedRHC, kUnweighted);
  //NoExtrapPredictionGenerator genNumu_FHC_IZZLE(axNumuEnergy, kIsNumuSelectedFHC, kUnweighted);
  //NoExtrapPredictionGenerator genNumu_RHC_IZZLE(axNumuEnergy, kIsNumuSelectedRHC, kUnweighted);

  PredictionInterp interpGenNue_FHC_IZZLE(systematicsVector, calc, genNue_FHC_IZZLE, loadersFHC, kNoShift, PredictionInterp::kSplitBySign);
  //PredictionInterp interpGenNue_RHC_IZZLE(systematicsVector, calc, *genNue_RHC_IZZLE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);
  //PredictionInterp interpGenNumu_FHC_IZZLE(systematicsVector, calc, *genNumu_FHC_IZZLE, loadersFHC, kNoShift, PredictionInterp::kSplitBySign);
  //PredictionInterp interpGenNumu_RHC_IZZLE(systematicsVector, calc, *genNumu_RHC_IZZLE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);

  std::cout << "Filling spectra..." << std::endl;
  loadersFHC.Go();
  //loadersRHC.Go();

  std::cout << "Saving to file: " << outputFileName << std::endl;
  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");

  interpGenNue_FHC_IZZLE.SaveTo(outputFile, "interpGenNue_FHC_IZZLE");
  //interpGenNue_RHC_IZZLE.SaveTo(outputFile, "interpGenNue_RHC_IZZLE");
  //interpGenNumu_FHC_IZZLE.SaveTo(outputFile, "interpGenNumu_FHC_IZZLE");
  //interpGenNumu_RHC_IZZLE.SaveTo(outputFile, "interpGenNumu_RHC_IZZLE");

  std::cout << "All done making state files..." << std::endl;
}

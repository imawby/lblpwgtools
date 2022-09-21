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
#include "CAFAna/Systs/DUNEFluxSysts.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

// also add in some RHC files!

using namespace ana;

const std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesFluxSystematicsSplitBySign.root";
//const std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesFluxSystematicsSplitBySign.root";
//const std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesEnergySystematicsSplitBySign.root";

void statisticalFluctuationStateFileProductionSystematics();

void statisticalFluctuationStateFileProductionSystematics()
{
   std::cout << "Reading caf files..." << std::endl;
   
  const std::string fnameNonSwapFHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/nu/caf_hist_nu.root";
  const std::string fnameNueSwapFHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/nue/caf_hist_nue.root";
  const std::string fnameTauSwapFHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/nutau/caf_hist_nutau.root";

  const std::string fnameNonSwapRHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/anu/caf_hist_anu.root";
  const std::string fnameNueSwapRHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/anue/caf_hist_anue.root";
  const std::string fnameTauSwapRHC = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/anutau/caf_hist_anutau.root";

  Loaders loadersFHC;
  loadersFHC.SetLoaderPath(fnameNonSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNonSwap);
  loadersFHC.SetLoaderPath(fnameNueSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNueSwap);
  loadersFHC.SetLoaderPath(fnameTauSwapFHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNuTauSwap);
  
  Loaders loadersRHC;
  loadersRHC.SetLoaderPath(fnameNonSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNonSwap);
  loadersRHC.SetLoaderPath(fnameNueSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNueSwap);
  loadersRHC.SetLoaderPath(fnameTauSwapRHC, caf::kFARDET, ana::Loaders::DataMC::kMC, ana::Loaders::SwappingConfig::kNuTauSwap);
  
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

  // Only flux systematics
  std::vector<ana::ISyst const *> systematicsVector = GetListOfSysts(true, false, false, false, true, false, true, false, 30, false);

  // Only xsec systematics
  //std::vector<ana::ISyst const *> systematicsVector = GetListOfSysts(false, true, false, false, true, false, true, false, 0, false);

  // Only energy systematics
  //std::vector<ana::ISyst const *> systematicsVector = GetListOfSysts(false, false, true, false, true, false, true, false, 0, false);

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
  NoExtrapPredictionGenerator genNue_RHC_IZZLE(axNueEnergy, kIsNueSelectedRHC, kUnweighted);
  NoExtrapPredictionGenerator genNumu_FHC_IZZLE(axNumuEnergy, kIsNumuSelectedFHC, kUnweighted);
  NoExtrapPredictionGenerator genNumu_RHC_IZZLE(axNumuEnergy, kIsNumuSelectedRHC, kUnweighted);

  const Cut kIsTrueNue([](const caf::SRProxy* sr)
  {
      return sr->isCC && abs(sr->nuPDGunosc) == 14 && abs(sr->nuPDG) == 12 && sr->nuPDG > 0;
  });

  const Cut kIsTrueAnue([](const caf::SRProxy* sr)
  {
      return sr->isCC && abs(sr->nuPDGunosc) == 14 && abs(sr->nuPDG) == 12 && sr->nuPDG < 0;
  });

  const Cut kIsTrueNumu([](const caf::SRProxy* sr)
  {
      return sr->isCC && abs(sr->nuPDGunosc) == 14 && abs(sr->nuPDG) == 14 && sr->nuPDG > 0;
  });

  const Cut kIsTrueAnumu([](const caf::SRProxy* sr)
  {
      return sr->isCC && abs(sr->nuPDGunosc) == 14 && abs(sr->nuPDG) == 14 && sr->nuPDG < 0;
  });

  NoExtrapPredictionGenerator genNue_FHC_TRUE(axNueEnergy, kIsTrueNue, kUnweighted);
  NoExtrapPredictionGenerator genNue_RHC_TRUE(axNueEnergy, kIsTrueAnue, kUnweighted);
  NoExtrapPredictionGenerator genNumu_FHC_TRUE(axNumuEnergy, kIsTrueNumu, kUnweighted);
  NoExtrapPredictionGenerator genNumu_RHC_TRUE(axNumuEnergy, kIsTrueAnumu, kUnweighted);

  NoExtrapPredictionGenerator genNue_FHC_CVN(axNueEnergy, kPassFD_CVN_NUE, kUnweighted);
  NoExtrapPredictionGenerator genNue_RHC_CVN(axNueEnergy, kPassFD_CVN_NUE, kUnweighted);
  NoExtrapPredictionGenerator genNumu_FHC_CVN(axNumuEnergy, kPassFD_CVN_NUMU, kUnweighted);
  NoExtrapPredictionGenerator genNumu_RHC_CVN(axNumuEnergy, kPassFD_CVN_NUMU, kUnweighted);

  PredictionInterp interpGenNue_FHC_IZZLE(systematicsVector, calc, genNue_FHC_IZZLE, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNue_RHC_IZZLE(systematicsVector, calc, genNue_RHC_IZZLE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);
  PredictionInterp interpGenNumu_FHC_IZZLE(systematicsVector, calc, genNumu_FHC_IZZLE, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNumu_RHC_IZZLE(systematicsVector, calc, genNumu_RHC_IZZLE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);

  PredictionInterp interpGenNue_FHC_TRUE(systematicsVector, calc, genNue_FHC_TRUE, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNue_RHC_TRUE(systematicsVector, calc, genNue_RHC_TRUE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);
  PredictionInterp interpGenNumu_FHC_TRUE(systematicsVector, calc, genNumu_FHC_TRUE, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNumu_RHC_TRUE(systematicsVector, calc, genNumu_RHC_TRUE, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);

  PredictionInterp interpGenNue_FHC_CVN(systematicsVector, calc, genNue_FHC_CVN, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNue_RHC_CVN(systematicsVector, calc, genNue_RHC_CVN, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);
  PredictionInterp interpGenNumu_FHC_CVN(systematicsVector, calc, genNumu_FHC_CVN, loadersFHC, kNoShift, PredictionInterp::kCombineSigns);
  PredictionInterp interpGenNumu_RHC_CVN(systematicsVector, calc, genNumu_RHC_CVN, loadersRHC, kNoShift, PredictionInterp::kSplitBySign);

  std::cout << "Filling spectra..." << std::endl;
  loadersFHC.Go();
  loadersRHC.Go();

  std::cout << "Saving to file: " << outputFileName << std::endl;
  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");

  interpGenNue_FHC_IZZLE.SaveTo(outputFile, "interpGenNue_FHC_IZZLE");
  interpGenNue_RHC_IZZLE.SaveTo(outputFile, "interpGenNue_RHC_IZZLE");
  interpGenNumu_FHC_IZZLE.SaveTo(outputFile, "interpGenNumu_FHC_IZZLE");
  interpGenNumu_RHC_IZZLE.SaveTo(outputFile, "interpGenNumu_RHC_IZZLE");

  interpGenNue_FHC_TRUE.SaveTo(outputFile, "interpGenNue_FHC_TRUE");
  interpGenNue_RHC_TRUE.SaveTo(outputFile, "interpGenNue_RHC_TRUE");
  interpGenNumu_FHC_TRUE.SaveTo(outputFile, "interpGenNumu_FHC_TRUE");
  interpGenNumu_RHC_TRUE.SaveTo(outputFile, "interpGenNumu_RHC_TRUE");

  interpGenNue_FHC_CVN.SaveTo(outputFile, "interpGenNue_FHC_CVN");
  interpGenNue_RHC_CVN.SaveTo(outputFile, "interpGenNue_RHC_CVN");
  interpGenNumu_FHC_CVN.SaveTo(outputFile, "interpGenNumu_FHC_CVN");
  interpGenNumu_RHC_CVN.SaveTo(outputFile, "interpGenNumu_RHC_CVN");

  std::cout << "All done making state files..." << std::endl;
}

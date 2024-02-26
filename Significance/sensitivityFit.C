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
#include "TTree.h"

// New includes for this macro
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Systs/AnaSysts.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/ISyst.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TRandom3.h"

#include <random>

using namespace ana;

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/realRecoCAF/trueCounts/StateFilesNoSystematicsSplitBySign_RecoEnergyRecoSelection.root";

//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/StateFilesNoSystematicsSplitBySign_RecoEnergyTrueSelection.root";
//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/StateFilesNoSystematicsSplitBySignTRUE.root";
//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/StateFilesNoSystematicsSplitBySign.root";
//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/StateFilesNoBinsNoSystematicsSplitBySign_STANDARD.root";
///storage/epp2/phrsnt/lblpwgtools/realRecoCAF/noSystematics/StateFilesNoSystematicsSplitBySign.root";

void sensitivityFit();

std::vector<Spectrum> Get_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);
void GetCount(const Spectrum &spectrum, const float pot, double &lowEnergyCount, double &middleEnergyCount, double &highEnergyCount, double &totalCount);
void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition, 
		double &bestChiSquared_nueFHC, double &bestChiSquared_numuFHC, double &bestChiSquared_nueRHC, double &bestChiSquared_numuRHC,
		double &bestChiSquared_FHC, double &bestChiSquared_RHC);
void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);

int N_TEST_DELTA_CP_VALUES = 200; // 200

void sensitivityFit()
{
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/TRUECountsTRUEEnergy_NO_SplitBySign.root";
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/TRUECountsRECOEnergy_NO_SplitBySign.root";
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/STANDARDCounts_NO_SplitBySign.root";
std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/realRecoCAF/trueCounts/STANDARDCounts_NO_SplitBySign.root";

  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/IZZLESensitivityPlotsNoBinsNoSystematicFit_NO_SplitBySign.root";

  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  
  PredictionInterp& interpGenNue_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_IZZLE").release();
  PredictionInterp& interpGenNue_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_IZZLE").release();
  PredictionInterp& interpGenNumu_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_IZZLE").release();
  PredictionInterp& interpGenNumu_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_IZZLE").release();
  
  //PredictionInterp& interpGenNue_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_TRUE").release();
  //PredictionInterp& interpGenNue_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_TRUE").release();
  //PredictionInterp& interpGenNumu_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_TRUE").release();
  //PredictionInterp& interpGenNumu_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_TRUE").release();

  //PredictionInterp& interpGenNue_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_CVN").release();
  //PredictionInterp& interpGenNue_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_CVN").release();
  //PredictionInterp& interpGenNumu_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_CVN").release();
  //PredictionInterp& interpGenNumu_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_CVN").release();

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_TRUE, &interpGenNue_RHC_TRUE, &interpGenNumu_FHC_TRUE, &interpGenNumu_RHC_TRUE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_CVN, &interpGenNue_RHC_CVN, &interpGenNumu_FHC_CVN, &interpGenNumu_RHC_CVN};

  inputFile->Close();

  const double pot = 3.5 * 1.1e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  double deltaCPValues;
  double chiSquaredCPCValues;
  double chiSquaredCPVValues;
  double chiSquaredValues;
  double chiSquaredNueFHC;
  double chiSquaredNumuFHC;
  double chiSquaredNueRHC;
  double chiSquaredNumuRHC;
  double chiSquaredFHC;
  double chiSquaredRHC;
  double bestFitDmSq32;
  double bestFitTh23;
  double bestFitTh12;
  double bestFitDmSq21;
  double bestFitTh13;
  double bestFitRho;
  double bestFitdCP;
  double nueLowEnergyCount, nueMiddleEnergyCount, nueHighEnergyCount, nueTotalCount;
  double anueLowEnergyCount, anueMiddleEnergyCount, anueHighEnergyCount, anueTotalCount;
  double numuLowEnergyCount, numuMiddleEnergyCount, numuHighEnergyCount, numuTotalCount;
  double anumuLowEnergyCount, anumuMiddleEnergyCount, anumuHighEnergyCount, anumuTotalCount;
  std::vector<double> nueTotalCountBin;
  std::vector<double> numuTotalCountBin;

  tree->Branch("deltaCPValues", &deltaCPValues);
  tree->Branch("chiSquaredValues" , &chiSquaredValues);
  tree->Branch("chiSquaredCPVValues" , &chiSquaredCPVValues);
  tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
  tree->Branch("chiSquaredNueFHC", &chiSquaredNueFHC);
  tree->Branch("chiSquaredNumuFHC", &chiSquaredNumuFHC);
  tree->Branch("chiSquaredNueRHC", &chiSquaredNueRHC);
  tree->Branch("chiSquaredNumuRHC", &chiSquaredNumuRHC);
  tree->Branch("chiSquaredFHC", &chiSquaredFHC);
  tree->Branch("chiSquaredRHC", &chiSquaredRHC);
  tree->Branch("bestFitDmSq32", &bestFitDmSq32);
  tree->Branch("bestFitTh23", &bestFitTh23);
  tree->Branch("bestFitTh12", &bestFitTh12);
  tree->Branch("bestFitDmSq21", &bestFitDmSq21);
  tree->Branch("bestFitTh13", &bestFitTh13);
  tree->Branch("bestFitRho", &bestFitRho);
  tree->Branch("bestFitdCP", &bestFitdCP);
  tree->Branch("nueLowEnergyCount", &nueLowEnergyCount);
  tree->Branch("nueMiddleEnergyCount", &nueMiddleEnergyCount);
  tree->Branch("nueHighEnergyCount", &nueHighEnergyCount);
  tree->Branch("nueTotalCount", &nueTotalCount);
  tree->Branch("anueLowEnergyCount", &anueLowEnergyCount);
  tree->Branch("anueMiddleEnergyCount", &anueMiddleEnergyCount);
  tree->Branch("anueHighEnergyCount", &anueHighEnergyCount);
  tree->Branch("anueTotalCount", &anueTotalCount);
  tree->Branch("numuLowEnergyCount", &numuLowEnergyCount);
  tree->Branch("numuMiddleEnergyCount", &numuMiddleEnergyCount);
  tree->Branch("numuHighEnergyCount", &numuHighEnergyCount);
  tree->Branch("numuTotalCount", &numuTotalCount);
  tree->Branch("anumuLowEnergyCount", &anumuLowEnergyCount);
  tree->Branch("anumuMiddleEnergyCount", &anumuMiddleEnergyCount);
  tree->Branch("anumuHighEnergyCount", &anumuHighEnergyCount);
  tree->Branch("anumuTotalCount", &anumuTotalCount);
  tree->Branch("nueTotalCountBin", &nueTotalCountBin);
  tree->Branch("numuTotalCountBin", &numuTotalCountBin);

  for (int j = 0; j < nTestCPValues; ++j)
  {
    nueTotalCountBin.clear();
    numuTotalCountBin.clear();

    std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;

    const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);
    deltaCPValues = trueDeltaCP;

    // Get the prediction
    std::vector<Spectrum> predictionVector;
    predictionVector = Get_NO(predictionGenerators, trueDeltaCP, pot);

    nueLowEnergyCount = 0; nueMiddleEnergyCount = 0; nueHighEnergyCount = 0; nueTotalCount = 0;
    anueLowEnergyCount = 0; anueMiddleEnergyCount = 0; anueHighEnergyCount = 0; anueTotalCount = 0;
    numuLowEnergyCount = 0; numuMiddleEnergyCount = 0; numuHighEnergyCount = 0; numuTotalCount = 0;
    anumuLowEnergyCount = 0; anumuMiddleEnergyCount = 0; anumuHighEnergyCount = 0; anumuTotalCount = 0;

    GetCount(predictionVector[0], pot, nueLowEnergyCount, nueMiddleEnergyCount, nueHighEnergyCount, nueTotalCount);
    GetCount(predictionVector[1], pot, anueLowEnergyCount, anueMiddleEnergyCount, anueHighEnergyCount, anueTotalCount);
    GetCount(predictionVector[2], pot, numuLowEnergyCount, numuMiddleEnergyCount, numuHighEnergyCount, numuTotalCount);
    GetCount(predictionVector[3], pot, anumuLowEnergyCount, anumuMiddleEnergyCount, anumuHighEnergyCount, anumuTotalCount);

    TH1 * histNue = predictionVector[0].ToTH1(pot);
    int nBinsNue = histNue->GetNbinsX();

    for (int i = 1; i <=nBinsNue; ++i)
    {
      const double binContent = histNue->GetBinContent(i);
      nueTotalCountBin.push_back(binContent);
    }

    TH1 * histNumu = predictionVector[3].ToTH1(pot);
    int nBinsNumu = histNumu->GetNbinsX();

    for (int i = 1; i <=nBinsNumu; ++i)
    {
      const double binContent = histNumu->GetBinContent(i);
      numuTotalCountBin.push_back(binContent);
    }

    // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
    double bestChiSquaredCPC(std::numeric_limits<float>::max());
    double bestChiSquaredCPV(std::numeric_limits<float>::max());
    double bestChiSquaredNueFHC(std::numeric_limits<float>::max());
    double bestChiSquaredNumuFHC(std::numeric_limits<float>::max());
    double bestChiSquaredNueRHC(std::numeric_limits<float>::max());
    double bestChiSquaredNumuRHC(std::numeric_limits<float>::max());
    double bestChiSquaredFHC(std::numeric_limits<float>::max());
    double bestChiSquaredRHC(std::numeric_limits<float>::max());

    std::map<std::string, float> bestFitPosition_CPC, bestFitPosition_CPV;

    // CPC Fits
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC, 
	       bestChiSquaredNueFHC, bestChiSquaredNumuFHC, bestChiSquaredNueRHC, bestChiSquaredNumuRHC, bestChiSquaredFHC, bestChiSquaredRHC);
    
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC,
	       bestChiSquaredNueFHC, bestChiSquaredNumuFHC, bestChiSquaredNueRHC, bestChiSquaredNumuRHC, bestChiSquaredFHC, bestChiSquaredRHC);

    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC,
	       bestChiSquaredNueFHC, bestChiSquaredNumuFHC, bestChiSquaredNueRHC, bestChiSquaredNumuRHC, bestChiSquaredFHC, bestChiSquaredRHC);

    /*
    std::cout << "Performing CPC fits..." << std::endl;
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

    std::cout << "Performing CPV fits..." << std::endl;
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, false, false, bestChiSquaredCPV, bestFitPosition_CPV);

    PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), true, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), true, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), false, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), false, false, false, bestChiSquaredCPV, bestFitPosition_CPV);

    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, false, false, bestChiSquaredCPV, bestFitPosition_CPV);

    PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), true, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), true, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), false, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), false, false, false, bestChiSquaredCPV, bestFitPosition_CPV);

    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, true, false, bestChiSquaredCPV, bestFitPosition_CPV);
    PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, false, false, bestChiSquaredCPV, bestFitPosition_CPV);
    */
    chiSquaredCPCValues = bestChiSquaredCPC;
    chiSquaredCPVValues = bestChiSquaredCPV;
    chiSquaredValues = bestChiSquaredCPC - bestChiSquaredCPV;
    chiSquaredNueFHC = bestChiSquaredNueFHC;
    chiSquaredNumuFHC = bestChiSquaredNumuFHC;
    chiSquaredNueRHC = bestChiSquaredNueRHC;
    chiSquaredNumuRHC = bestChiSquaredNumuRHC;
    chiSquaredFHC = bestChiSquaredFHC;
    chiSquaredRHC = bestChiSquaredRHC;
    bestFitDmSq32 = bestFitPosition_CPC.find("kFitDmSq32Scaled") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq32Scaled") : -999.0;
    bestFitTh23 = bestFitPosition_CPC.find("kFitSinSqTheta23") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSqTheta23") : -999.0;
    bestFitTh12 = bestFitPosition_CPC.find("kFitSinSq2Theta12") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSq2Theta12") : -999.0;
    bestFitDmSq21 = bestFitPosition_CPC.find("kFitDmSq21") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq21") : -999.0;
    bestFitTh13 = bestFitPosition_CPC.find("kFitTheta13") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitTheta13") : -999.0;
    bestFitRho = bestFitPosition_CPC.find("kFitRho") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitRho") : -999.0;
    bestFitdCP = bestFitPosition_CPC.find("kFitDeltaInPiUnits") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDeltaInPiUnits") : -999.0;

    tree->Fill();
  }

  outputFile->WriteObject(tree, "tree");    
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot)
{
    // Make oscillation calc
    osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
    calc->SetdCP(deltaCP);

    // Create experiments
    std::vector<Spectrum> predictionVector;

    for (const PredictionInterp *const predictionGenerator : predictionGenerators)
    {
        const Spectrum prediction(predictionGenerator->Predict(calc).AsimovData(pot));
        predictionVector.push_back(prediction);
    }

    return predictionVector;
}

////////////////////////////////////////////////////////////////////////////////////////////

void GetCount(const Spectrum &spectrum, const float pot, double &lowEnergyCount, double &middleEnergyCount, double &highEnergyCount, double &totalCount)
{
  TH1 * hist = spectrum.ToTH1(pot);
  int nBins = hist->GetNbinsX();

  for (int i = 1; i <=nBins; ++i)
  {
    const double binCenter = hist->GetBinCenter(i);
    const double binContent = hist->GetBinContent(i);

    if (binCenter < 1.50)
      lowEnergyCount += binContent;
    else if ((binCenter > 1.50) && (binCenter < 3.50))
      middleEnergyCount += binContent;
    else if ((binCenter > 3.50) && (binCenter < 6.0))
      highEnergyCount += binContent;

    totalCount += binContent;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition, 
		double &bestChiSquared_nueFHC, double &bestChiSquared_numuFHC, double &bestChiSquared_nueRHC, double &bestChiSquared_numuRHC,
		double &bestChiSquared_FHC, double &bestChiSquared_RHC)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  //osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
  calc->SetdCP(deltaCPSeed);

  // Get fit variables
  std::vector<const IFitVar*> fitVariables;/* =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};*/
					   
  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});
  MultiExperiment multiExperiment_nueFHC({&experimentNueFHC});
  MultiExperiment multiExperiment_numuFHC({&experimentNumuFHC});
  MultiExperiment multiExperiment_nueRHC({&experimentNueRHC});
  MultiExperiment multiExperiment_numuRHC({&experimentNumuRHC});
  MultiExperiment multiExperiment_FHC({&experimentNueFHC, &experimentNumuFHC});
  MultiExperiment multiExperiment_RHC({&experimentNueRHC, &experimentNumuRHC});

  MinuitFitter fit(&multiExperiment, fitVariables);
  float chiSquared = fit.Fit(calc)->EvalMetricVal();

   if (chiSquared < bestChiSquared)
   {
       bestFitPosition["kFitDmSq32Scaled"] = calc->GetDmsq32();
       bestFitPosition["kFitSinSqTheta23"] = calc->GetTh23();
       bestFitPosition["kFitSinSq2Theta12"] = calc->GetTh12();
       bestFitPosition["kFitDmSq21"] = calc->GetDmsq21();
       bestFitPosition["kFitTheta13"] = calc->GetTh13();
       bestFitPosition["kFitRho"] = calc->GetRho();
       bestFitPosition["kFitDeltaInPiUnits"] = calc->GetdCP();

       bestChiSquared = chiSquared;
   }

  MinuitFitter fit_nueFHC(&multiExperiment_nueFHC, fitVariables);
  float chiSquared_nueFHC = fit_nueFHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_nueFHC < bestChiSquared_nueFHC)
    bestChiSquared_nueFHC = chiSquared_nueFHC;

  MinuitFitter fit_numuFHC(&multiExperiment_numuFHC, fitVariables);
  float chiSquared_numuFHC = fit_numuFHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_numuFHC < bestChiSquared_numuFHC)
    bestChiSquared_numuFHC = chiSquared_numuFHC;

  MinuitFitter fit_nueRHC(&multiExperiment_nueRHC, fitVariables);
  float chiSquared_nueRHC = fit_nueRHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_nueRHC < bestChiSquared_nueRHC)
    bestChiSquared_nueRHC = chiSquared_nueRHC;

  MinuitFitter fit_numuRHC(&multiExperiment_numuRHC, fitVariables);
  float chiSquared_numuRHC = fit_numuRHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_numuRHC < bestChiSquared_numuRHC)
    bestChiSquared_numuRHC = chiSquared_numuRHC;

  MinuitFitter fit_FHC(&multiExperiment_FHC, fitVariables);
  float chiSquared_FHC = fit_FHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_FHC < bestChiSquared_FHC)
    bestChiSquared_FHC = chiSquared_FHC;

  MinuitFitter fit_RHC(&multiExperiment_RHC, fitVariables);
  float chiSquared_RHC = fit_RHC.Fit(calc)->EvalMetricVal();

  if (chiSquared_RHC < bestChiSquared_RHC)
    bestChiSquared_RHC = chiSquared_RHC;
}

////////////////////////////////////////////////////////////////////////////////////////////

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
  const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  //osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
  calc->SetdCP(deltaCPSeed);

  // Get fit variables
  std::vector<const IFitVar*> fitVariables =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};
					   
  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});

  MinuitFitter fit(&multiExperiment, fitVariables);
  float chiSquared = fit.Fit(calc)->EvalMetricVal();

   if (chiSquared < bestChiSquared)
   {
       bestFitPosition["kFitDmSq32Scaled"] = calc->GetDmsq32();
       bestFitPosition["kFitSinSqTheta23"] = calc->GetTh23();
       bestFitPosition["kFitSinSq2Theta12"] = calc->GetTh12();
       bestFitPosition["kFitDmSq21"] = calc->GetDmsq21();
       bestFitPosition["kFitTheta13"] = calc->GetTh13();
       bestFitPosition["kFitRho"] = calc->GetRho();
       bestFitPosition["kFitDeltaInPiUnits"] = calc->GetdCP();

       bestChiSquared = chiSquared;
   }
}

////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////

  /*std::cout << "after: " << std::endl;
  std::cout << "rho: " << calc->GetRho() << std::endl;
  std::cout << "dmsq21: " << calc->GetDmsq21() << std::endl;
  std::cout << "dmsq32: " << calc->GetDmsq32() << std::endl;
  std::cout << "theta12: " << calc->GetTh12() << std::endl;
  std::cout << "theta13: " << calc->GetTh13() << std::endl;
  std::cout << "theta23: " << calc->GetTh23() << std::endl;
  std::cout << "DeltaCP: " << calc->GetdCP() << std::endl;*/

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

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/StateFilesNoSystematicsSplitBySignTRUE.root";

void sensitivityFit();

std::vector<Spectrum> Get_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);
void GetCount(const Spectrum &spectrum, const float pot, double &lowEnergyCount, double &middleEnergyCount, double &highEnergyCount, double &totalCount);
void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);

int N_TEST_DELTA_CP_VALUES = 200; // 200

void sensitivityFit()
{
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/trueCounts/TRUECounts_NO_SplitBySign.root";

  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  //PredictionInterp& interpGenNue_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_IZZLE").release();
  //PredictionInterp& interpGenNue_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_IZZLE").release();
  //PredictionInterp& interpGenNumu_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_IZZLE").release();
  //PredictionInterp& interpGenNumu_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_IZZLE").release();

  PredictionInterp& interpGenNue_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_TRUE").release();
  PredictionInterp& interpGenNue_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_TRUE").release();
  PredictionInterp& interpGenNumu_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_TRUE").release();
  PredictionInterp& interpGenNumu_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_TRUE").release();

  //PredictionInterp& interpGenNue_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_CVN").release();
  //PredictionInterp& interpGenNue_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_CVN").release();
  //PredictionInterp& interpGenNumu_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_CVN").release();
  //PredictionInterp& interpGenNumu_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_CVN").release();

  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};
  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_TRUE, &interpGenNue_RHC_TRUE, &interpGenNumu_FHC_TRUE, &interpGenNumu_RHC_TRUE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_CVN, &interpGenNue_RHC_CVN, &interpGenNumu_FHC_CVN, &interpGenNumu_RHC_CVN};

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  double deltaCPValues;
  double chiSquaredCPCValues;
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

  tree->Branch("deltaCPValues", &deltaCPValues);
  tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
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


  for (int j = 0; j < nTestCPValues; ++j)
  {
    nueTotalCountBin.clear();

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

  TH1 * hist = predictionVector[0].ToTH1(pot);
  int nBins = hist->GetNbinsX();

  for (int i = 1; i <=nBins; ++i)
  {
    const double binContent = hist->GetBinContent(i);
    nueTotalCountBin.push_back(binContent);
  }


    // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
    double bestChiSquaredCPC(std::numeric_limits<float>::max());
    std::map<std::string, float> bestFitPosition_CPC;

    // CPC Fits
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

    chiSquaredCPCValues = bestChiSquaredCPC;
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
    calc->SetTh23(40.3*3.14/180);
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
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  //osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
  calc->SetdCP(deltaCPSeed);
  calc->SetTh23(40.3*3.14/180);

  // Get fit variables
  std::vector<const IFitVar*> fitVariables; /*=
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};
					    */
  
  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  //const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  //const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  //const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  //MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});
  MultiExperiment multiExperiment({&experimentNueFHC});
  //MultiExperiment multiExperiment({&experimentNueFHC, &experimentNumuFHC});

  float chiSquared(std::numeric_limits<float>::max());

  MinuitFitter fit(&multiExperiment, fitVariables);
  chiSquared = fit.Fit(calc)->EvalMetricVal();

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

  /*std::cout << "after: " << std::endl;
  std::cout << "rho: " << calc->GetRho() << std::endl;
  std::cout << "dmsq21: " << calc->GetDmsq21() << std::endl;
  std::cout << "dmsq32: " << calc->GetDmsq32() << std::endl;
  std::cout << "theta12: " << calc->GetTh12() << std::endl;
  std::cout << "theta13: " << calc->GetTh13() << std::endl;
  std::cout << "theta23: " << calc->GetTh23() << std::endl;
  std::cout << "DeltaCP: " << calc->GetdCP() << std::endl;*/

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

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesEnergySystematicsSplitBySign.root";

void sensitivityFitTh23();

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double trueTh23, const double deltaCPSeed, 
  const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);
std::vector<Spectrum> Get_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double trueTh23, const double deltaCP, const float pot);

int N_TEST_DELTA_CP_VALUES = 200; // 200

void sensitivityFitTh23()
{
  // For th23 studies
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/th23/IZZLESensitivityPlotsNoSystematics_NO_SplitBySign_FHC.root";

  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESpectraPlotsSystematics_NO_SplitBySign_LOW.root";

  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  PredictionInterp& interpGenNue_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_IZZLE").release();
  PredictionInterp& interpGenNue_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_IZZLE").release();
  PredictionInterp& interpGenNumu_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_IZZLE").release();
  PredictionInterp& interpGenNumu_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_IZZLE").release();

  //PredictionInterp& interpGenNue_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_CVN").release();
  //PredictionInterp& interpGenNue_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_CVN").release();
  //PredictionInterp& interpGenNumu_FHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_CVN").release();
  //PredictionInterp& interpGenNumu_RHC_CVN = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_CVN").release();

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_CVN, &interpGenNue_RHC_CVN, &interpGenNumu_FHC_CVN, &interpGenNumu_RHC_CVN};

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  double th23Values;
  std::vector<double> chiSquaredCPCValues;
  std::vector<double> bestFitDmSq32;
  std::vector<double> bestFitTh23;
  std::vector<double> bestFitTh12;
  std::vector<double> bestFitDmSq21;
  std::vector<double> bestFitTh13;
  std::vector<double> bestFitRho;
  std::vector<double> bestFitdCP;

  tree->Branch("th23Values", &th23Values);
  tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
  tree->Branch("bestFitDmSq32", &bestFitDmSq32);
  tree->Branch("bestFitTh23", &bestFitTh23);
  tree->Branch("bestFitTh12", &bestFitTh12);
  tree->Branch("bestFitDmSq21", &bestFitDmSq21);
  tree->Branch("bestFitTh13", &bestFitTh13);
  tree->Branch("bestFitRho", &bestFitRho);
  tree->Branch("bestFitdCP", &bestFitdCP);

  int nTh23Values = 100;

  for (int i = 0; i < nTh23Values; ++i)
  {
    // th23
    double stepSizeDegrees = (180 - 0) / nTh23Values;
    double th23Degrees = 0 + (i * stepSizeDegrees);

    // dm31
    //double stepSizeDegrees = (2.625 - 2.427) / nTh23Values;
    //double th23Degrees = 2.427 + (i * stepSizeDegrees);

    th23Values = th23Degrees;

    chiSquaredCPCValues.clear();
    bestFitDmSq32.clear();
    bestFitTh23.clear();
    bestFitTh12.clear();
    bestFitDmSq21.clear();
    bestFitTh13.clear();
    bestFitRho.clear();
    bestFitdCP.clear();

    for (int j = 0; j < nTestCPValues; ++j)
    {
      std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;

      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);
      const double trueTh23 = th23Degrees * TMath::Pi() / 180.0; 
      //const double trueTh23 = th23Degrees * 0.001; 

      // Get the prediction
      std::vector<Spectrum> predictionVector;
      predictionVector = Get_NO(predictionGenerators, trueTh23, trueDeltaCP, pot);

      // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
      double bestChiSquaredCPC(std::numeric_limits<float>::max());
      std::map<std::string, float> bestFitPosition_CPC;

      // CPC Fits
      std::cout << "Performing CPC fits..." << std::endl;
      PerformFit(predictionGenerators, predictionVector, trueTh23, 0.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 0.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 0.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 0.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

      PerformFit(predictionGenerators, predictionVector, trueTh23, 1.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 1.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 1.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 1.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

      PerformFit(predictionGenerators, predictionVector, trueTh23, 2.0 * TMath::Pi(), true, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 2.0 * TMath::Pi(), true, false, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 2.0 * TMath::Pi(), false, true, true, bestChiSquaredCPC, bestFitPosition_CPC);
      PerformFit(predictionGenerators, predictionVector, trueTh23, 2.0 * TMath::Pi(), false, false, true, bestChiSquaredCPC, bestFitPosition_CPC);

      chiSquaredCPCValues.push_back(bestChiSquaredCPC);

      bestFitDmSq32.push_back(bestFitPosition_CPC.find("kFitDmSq32Scaled") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq32Scaled") : -999.0);
      bestFitTh23.push_back(bestFitPosition_CPC.find("kFitSinSqTheta23") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSqTheta23") : -999.0);
      bestFitTh12.push_back(bestFitPosition_CPC.find("kFitSinSq2Theta12") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSq2Theta12") : -999.0);
      bestFitDmSq21.push_back(bestFitPosition_CPC.find("kFitDmSq21") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq21") : -999.0);
      bestFitTh13.push_back(bestFitPosition_CPC.find("kFitTheta13") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitTheta13") : -999.0);
      bestFitRho.push_back(bestFitPosition_CPC.find("kFitRho") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitRho") : -999.0);
      bestFitdCP.push_back(trueDeltaCP);
  }

      tree->Fill();
  }

  outputFile->WriteObject(tree, "tree");    
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double trueTh23, const double deltaCP, const float pot)
{
    // Make oscillation calc
    osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
    calc->SetdCP(deltaCP);

    // for th23 investigation
    calc->SetTh23(trueTh23);
    //calc->SetDmsq32(trueTh23); 

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

 void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double trueTh23, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  calc->SetdCP(deltaCPSeed);

  calc->SetTh23(trueTh23);
  //calc->SetDmsq32(trueTh23);

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
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  //MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});
  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNumuFHC});

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

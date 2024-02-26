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
#include "CAFAna/Analysis/common_fit_definitions.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TRandom3.h"

#include <random>

using namespace ana;

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesEnergySystematicsSplitBySign.root";

void sensitivityFitNoSystematics();

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);
std::vector<Spectrum> Get_NoSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot, std::map<std::string, float> &truePosition);

std::default_random_engine generator;

int N_TEST_DELTA_CP_VALUES = 200;
int N_THROWS = 500;

void sensitivityFitNoSystematics()
{
  // Change fit to throws
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/IZZLESensitivityPlotsNoSystematicThrow_NO_SplitBySign.root";
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/CVNSensitivityPlotsNoSystematicFit_NO_SplitBySign.root";

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
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNumu_FHC_IZZLE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_CVN, &interpGenNue_RHC_CVN, &interpGenNumu_FHC_CVN, &interpGenNumu_RHC_CVN};

  inputFile->Close();

  
  const double pot = 3.5 * 1.1e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  unsigned int nThrows(N_THROWS);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  double deltaCPValues;
  std::vector<double> chiSquaredCPVValues;
  std::vector<double> bestFitDmSq32_CPV;
  std::vector<double> bestFitTh23_CPV;
  std::vector<double> bestFitTh12_CPV;
  std::vector<double> bestFitDmSq21_CPV;
  std::vector<double> bestFitTh13_CPV;
  std::vector<double> bestFitRho_CPV;
  std::vector<double> bestFitdCP_CPV;

  std::vector<double> trueDmSq32;
  std::vector<double> trueTh23;
  std::vector<double> trueTh12;
  std::vector<double> trueDmSq21;
  std::vector<double> trueTh13;
  std::vector<double> trueRho;
  std::vector<double> truedCP;

  tree->Branch("deltaCPValues", &deltaCPValues);
  tree->Branch("chiSquaredCPVValues", &chiSquaredCPVValues);
  tree->Branch("bestFitDmSq32_CPV", &bestFitDmSq32_CPV);
  tree->Branch("bestFitTh23_CPV", &bestFitTh23_CPV);
  tree->Branch("bestFitTh12_CPV", &bestFitTh12_CPV);
  tree->Branch("bestFitDmSq21_CPV", &bestFitDmSq21_CPV);
  tree->Branch("bestFitTh13_CPV", &bestFitTh13_CPV);
  tree->Branch("bestFitRho_CPV", &bestFitRho_CPV);
  tree->Branch("bestFitdCP_CPV", &bestFitdCP_CPV);

  tree->Branch("trueDmSq32", &trueDmSq32);
  tree->Branch("trueTh23", &trueTh23);
  tree->Branch("trueTh12", &trueTh12);
  tree->Branch("trueDmSq21", &trueDmSq21);
  tree->Branch("trueTh13", &trueTh13);
  tree->Branch("trueRho", &trueRho);
  tree->Branch("truedCP", &truedCP);

  for (int j = 0; j < nTestCPValues; ++j)
  {
      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);

      deltaCPValues = trueDeltaCP;

      chiSquaredCPVValues.clear();
      bestFitDmSq32_CPV.clear();
      bestFitTh23_CPV.clear();
      bestFitTh12_CPV.clear();
      bestFitDmSq21_CPV.clear();
      bestFitTh13_CPV.clear();
      bestFitRho_CPV.clear();
      bestFitdCP_CPV.clear();

      trueDmSq32.clear();
      trueTh23.clear();
      trueTh12.clear();
      trueDmSq21.clear();
      trueTh13.clear();
      trueRho.clear();
      truedCP.clear();

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i + 1) << "/" << std::to_string(nThrows) << std::endl;

          // Get the prediction
          std::vector<Spectrum> predictionVector;
	  std::map<std::string, float> truePosition;

	  predictionVector = Get_NoSys_NO(predictionGenerators, trueDeltaCP, pot, truePosition);

	  trueDmSq32.push_back(truePosition.find("kFitDmSq32Scaled") != truePosition.end() ? truePosition.at("kFitDmSq32Scaled") : -999.0);
	  trueTh23.push_back(truePosition.find("kFitSinSqTheta23") != truePosition.end() ? truePosition.at("kFitSinSqTheta23") : -999.0);
	  trueTh12.push_back(truePosition.find("kFitSinSq2Theta12") != truePosition.end() ? truePosition.at("kFitSinSq2Theta12") : -999.0);
	  trueDmSq21.push_back(truePosition.find("kFitDmSq21") != truePosition.end() ? truePosition.at("kFitDmSq21") : -999.0);
	  trueTh13.push_back(truePosition.find("kFitTheta13") != truePosition.end() ? truePosition.at("kFitTheta13") : -999.0);
	  trueRho.push_back(truePosition.find("kFitRho") != truePosition.end() ? truePosition.at("kFitRho") : -999.0);
	  truedCP.push_back(truePosition.find("kFitDeltaInPiUnits") != truePosition.end() ? truePosition.at("kFitDeltaInPiUnits") : -999.0);

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

	  std::map<std::string, float> bestFitPosition_CPV;

	  PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), true, false, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), false, false, bestChiSquaredCPV, bestFitPosition_CPV);

	  PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), true, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), true, false, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), false, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), false, false, bestChiSquaredCPV, bestFitPosition_CPV);

	  PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), true, false, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), false, false, bestChiSquaredCPV, bestFitPosition_CPV);

	  PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), true, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), true, false, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), false, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), false, false, bestChiSquaredCPV, bestFitPosition_CPV);

	  PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), true, false, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, true, bestChiSquaredCPV, bestFitPosition_CPV);
	  PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), false, false, bestChiSquaredCPV, bestFitPosition_CPV);

	  chiSquaredCPVValues.push_back(bestChiSquaredCPV);
	  bestFitDmSq32_CPV.push_back(bestFitPosition_CPV.find("kFitDmSq32Scaled") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDmSq32Scaled") : -999.0);
	  bestFitTh23_CPV.push_back(bestFitPosition_CPV.find("kFitSinSqTheta23") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitSinSqTheta23") : -999.0);
	  bestFitTh12_CPV.push_back(bestFitPosition_CPV.find("kFitSinSq2Theta12") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitSinSq2Theta12") : -999.0);
	  bestFitDmSq21_CPV.push_back(bestFitPosition_CPV.find("kFitDmSq21") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDmSq21") : -999.0);
	  bestFitTh13_CPV.push_back(bestFitPosition_CPV.find("kFitTheta13") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitTheta13") : -999.0);
	  bestFitRho_CPV.push_back(bestFitPosition_CPV.find("kFitRho") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitRho") : -999.0);
	  bestFitdCP_CPV.push_back(bestFitPosition_CPV.find("kFitDeltaInPiUnits") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDeltaInPiUnits") : -999.0);
      }  

      tree->Fill();
  }

  outputFile->WriteObject(tree, "tree");    
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_NoSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot, std::map<std::string, float> &truePosition)
{
    // Make oscillation calc throw (1 is for NO)
    std::vector<const IFitVar *> oscVars = GetOscVars("dmsq32:th23:th13", 1);

    osc::IOscCalcAdjustable* calc = ThrownWideOscCalc(1, oscVars, false);
    calc->SetdCP(deltaCP);

    truePosition["kFitDmSq32Scaled"] = calc->GetDmsq32();
    truePosition["kFitSinSqTheta23"] = calc->GetTh23();
    truePosition["kFitSinSq2Theta12"] = calc->GetTh12();
    truePosition["kFitDmSq21"] = calc->GetDmsq21();
    truePosition["kFitTheta13"] = calc->GetTh13();
    truePosition["kFitRho"] = calc->GetRho();
    truePosition["kFitDeltaInPiUnits"] = calc->GetdCP();

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

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  calc->SetdCP(deltaCPSeed);

  // Get fit variables
  std::vector<const IFitVar*> fitVariables =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho, &kFitDeltaInPiUnits};

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  /*
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);
  */

  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[1], predictionVector[1]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  //MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC});//, &penalty});
  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNumuFHC});//, &penalty});

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

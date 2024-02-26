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

//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/realRecoCAF/fullEstimate/StateFilesAllSystematicsSplitBySign_REAL.root";
const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesAllSystematicsSplitBySign_STANDARD.root";
///storage/epp2/phrsnt/lblpwgtools/realRecoStandardCAF/fullEstimate/StateFilesAllSystematicsSplitBySign.root";
///storage/epp2/phrsnt/lblpwgtools/standardCAF/fullEstimate/StateFilesAllSystematicsSplitBySign.root";
///storage/epp2/phrsnt/lblpwgtools/standardCAF/fullEstimate/StateFilesAllSystematicsSplitBySign.root";
//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/realRecoCAF/noSystematics/StateFilesNoSystematicsSplitBySign.root";

void sensitivityFitNoSystematicsOscFit();

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
  const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);
std::vector<Spectrum> Get_NoSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);

int N_TEST_DELTA_CP_VALUES = 200; // 200
bool MAKE_THROWS = false;
int N_THROWS = MAKE_THROWS ? 500 : 1; //

void sensitivityFitNoSystematicsOscFit()
{
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/realRecoCAF/noSystematics/IZZLESensitivityPlotsNoSystematicFit_NO_SplitBySign.root";
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/IZZLESensitivityPlotsNoSystematicFit_DM32_HalfPOT_NO_SplitBySign.root";
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/noSystematics/IZZLESensitivityPlotsNoSystematicFit_AllOsc_NO_SplitBySign.root";

  DUNEFluxSystVector fluxSystematics = GetDUNEFluxSysts(30, true, false);

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

  //PredictionInterp& interpGenNue_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_TRUE").release();
  //PredictionInterp& interpGenNue_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_TRUE").release();
  //PredictionInterp& interpGenNumu_FHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_TRUE").release();
  //PredictionInterp& interpGenNumu_RHC_TRUE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_TRUE").release();

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};
  //std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_CVN, &interpGenNue_RHC_CVN, &interpGenNumu_FHC_CVN, &interpGenNumu_RHC_CVN};
  //std::vector<const PredictionInterp*> predictionGenerators_TRUE = {&interpGenNue_FHC_TRUE, &interpGenNue_RHC_TRUE, &interpGenNumu_FHC_TRUE, &interpGenNumu_RHC_TRUE};

  inputFile->Close();

  //const double pot = 3.5 * 1.47e21 * 40/1.13;
  const double pot = 3.5 * 1.1e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  unsigned int nThrows(N_THROWS);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  double deltaCPValues;
  std::vector<double> chiSquaredCPCValues;
  std::vector<double> chiSquaredCPVValues;
  std::vector<double> chiSquaredValues;

  std::vector<double> bestFitDmSq32_CPC;
  std::vector<double> bestFitTh23_CPC;
  std::vector<double> bestFitTh12_CPC;
  std::vector<double> bestFitDmSq21_CPC;
  std::vector<double> bestFitTh13_CPC;
  std::vector<double> bestFitRho_CPC;
  std::vector<double> bestFitdCP_CPC;

  std::vector<double> bestFitDmSq32_CPV;
  std::vector<double> bestFitTh23_CPV;
  std::vector<double> bestFitTh12_CPV;
  std::vector<double> bestFitDmSq21_CPV;
  std::vector<double> bestFitTh13_CPV;
  std::vector<double> bestFitRho_CPV;
  std::vector<double> bestFitdCP_CPV;

  tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
  tree->Branch("chiSquaredCPVValues", &chiSquaredCPVValues);
  tree->Branch("chiSquaredValues", &chiSquaredValues);
  tree->Branch("deltaCPValues", &deltaCPValues);

  tree->Branch("bestFitDmSq32_CPC", &bestFitDmSq32_CPC);
  tree->Branch("bestFitTh23_CPC", &bestFitTh23_CPC);
  tree->Branch("bestFitTh12_CPC", &bestFitTh12_CPC);
  tree->Branch("bestFitDmSq21_CPC", &bestFitDmSq21_CPC);
  tree->Branch("bestFitTh13_CPC", &bestFitTh13_CPC);
  tree->Branch("bestFitRho_CPC", &bestFitRho_CPC);
  tree->Branch("bestFitdCP_CPC", &bestFitdCP_CPC);

  tree->Branch("bestFitDmSq32_CPV", &bestFitDmSq32_CPV);
  tree->Branch("bestFitTh23_CPV", &bestFitTh23_CPV);
  tree->Branch("bestFitTh12_CPV", &bestFitTh12_CPV);
  tree->Branch("bestFitDmSq21_CPV", &bestFitDmSq21_CPV);
  tree->Branch("bestFitTh13_CPV", &bestFitTh13_CPV);
  tree->Branch("bestFitRho_CPV", &bestFitRho_CPV);
  tree->Branch("bestFitdCP_CPV", &bestFitdCP_CPV);

  float highestSensitivity = 0;

  for (int j = 0; j < nTestCPValues; ++j)
  {
      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);

      //const double trueDeltaCP(1.5 * 3.14);

      deltaCPValues = trueDeltaCP;

      chiSquaredCPCValues.clear();
      chiSquaredCPVValues.clear();
      chiSquaredValues.clear();

      bestFitDmSq32_CPC.clear();
      bestFitTh23_CPC.clear();
      bestFitTh12_CPC.clear();
      bestFitDmSq21_CPC.clear();
      bestFitTh13_CPC.clear();
      bestFitRho_CPC.clear();
      bestFitdCP_CPC.clear();

      bestFitDmSq32_CPV.clear();
      bestFitTh23_CPV.clear();
      bestFitTh12_CPV.clear();
      bestFitDmSq21_CPV.clear();
      bestFitTh13_CPV.clear();
      bestFitRho_CPV.clear();
      bestFitdCP_CPV.clear();

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i + 1) << "/" << std::to_string(nThrows) << std::endl;

          // Get the prediction
          std::vector<Spectrum> predictionVector;

	  predictionVector = Get_NoSys_NO(predictionGenerators, trueDeltaCP, pot);

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          double bestChiSquaredCPC(std::numeric_limits<float>::max());
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

	  std::map<std::string, float> bestFitPosition_CPC;
	  std::map<std::string, float> bestFitPosition_CPV;

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

            // CPV Fits

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

          const float chiSquaredDifference(bestChiSquaredCPC - bestChiSquaredCPV);

	  if (chiSquaredDifference > highestSensitivity)
	    highestSensitivity = chiSquaredDifference;

          //std::cout << "////////////////////////////" << std::endl;
          //std::cout << "trueDeltaCP: " << (trueDeltaCP * 180 / 3.14) << std::endl;
          //std::cout << "bestChiSquaredCPC: " << bestChiSquaredCPC << std::endl;
          //std::cout << "bestChiSquaredCPV: " << bestChiSquaredCPV << std::endl;
          //std::cout << "chiSquaredDifference: " << chiSquaredDifference << std::endl;
          //std::cout << "significance: " << TMath::Sqrt(std::fabs(chiSquaredDifference)) << std::endl;
          //std::cout << "////////////////////////////" << std::endl;

	  chiSquaredCPCValues.push_back(bestChiSquaredCPC);
	  chiSquaredCPVValues.push_back(bestChiSquaredCPV);
          chiSquaredValues.push_back(chiSquaredDifference);

	  bestFitDmSq32_CPC.push_back(bestFitPosition_CPC.find("kFitDmSq32Scaled") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq32Scaled") : -999.0);
	  bestFitTh23_CPC.push_back(bestFitPosition_CPC.find("kFitSinSqTheta23") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSqTheta23") : -999.0);
	  bestFitTh12_CPC.push_back(bestFitPosition_CPC.find("kFitSinSq2Theta12") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSq2Theta12") : -999.0);
	  bestFitDmSq21_CPC.push_back(bestFitPosition_CPC.find("kFitDmSq21") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq21") : -999.0);
	  bestFitTh13_CPC.push_back(bestFitPosition_CPC.find("kFitTheta13") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitTheta13") : -999.0);
	  bestFitRho_CPC.push_back(bestFitPosition_CPC.find("kFitRho") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitRho") : -999.0);
	  bestFitdCP_CPC.push_back(bestFitPosition_CPC.find("kFitDeltaInPiUnits") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDeltaInPiUnits") : -999.0);

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

  std::cout << "highestSensitivity: " << highestSensitivity << std::endl;

  outputFile->WriteObject(tree, "tree");    
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_NoSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot)
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

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
		const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Set the oscillation calculator seeed
  int hie = (isPositiveHierarchy ? 1 : -1);
  int oct = (isHigherOctant ? 1 : -1);

  osc::IOscCalcAdjustable* calc = NuFitOscCalc(hie, oct);
  calc->SetdCP(deltaCPSeed);

  // Get fit variables

  /*
  std::vector<const IFitVar*> fitVariables =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};
  */

  /*
  std::vector<const IFitVar*> fitVariables =
    {&kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};
  */

  std::vector<const IFitVar*> fitVariables;


  //std::vector<const IFitVar*> fitVariables =
  //{&kFitDmSq32Scaled};

  //std::vector<const IFitVar*> fitVariables =
  //{&kFitSinSq2Theta12};

  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);
  //Penalizer_GlbLike penalty(hie, oct, false, false, false, 0);

  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});
  //MultiExperiment multiExperiment({&experimentNueFHC, &experimentNumuFHC, &penalty});

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

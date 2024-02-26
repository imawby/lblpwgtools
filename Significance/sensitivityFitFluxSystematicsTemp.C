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

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesFluxSystematicsSplitBySign.root";

void sensitivityFitFluxSystematicsTemp(float sigmaShift);

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
  const bool isHigherOctant, const bool isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);
std::vector<Spectrum> Get_Thrown_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);
std::vector<Spectrum> Get_FluxSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, std::vector<std::string> &systematicsToShift, const float sigmaShift, const double deltaCP, const float pot);
Spectrum Get_CentralValue_NO(const PredictionInterp &predictionGenerator, const double deltaCP, const float pot);
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

std::default_random_engine generator;

int N_TEST_DELTA_CP_VALUES = 200; // 10
bool MAKE_THROWS = false;
int N_THROWS = MAKE_THROWS ? 200 : 1; // 500
bool FIT_SYSTEMATICS = true;

void sensitivityFitFluxSystematicsTemp(float sigmaShift)
{
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsFluxSystematics" + std::to_string(sigmaShift) + "SigmaShift_NO_SplitBySign.root";
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotFittedFluxSystematics_NO_SplitBySign.root";
  //std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESpectraPlotsFluxSystematics_NO_SplitBySign_ZERO.root";

  DUNEFluxSystVector systematicsVector = GetDUNEFluxSysts(30, true, false);

  std::vector<std::string> systematicsToShift;

  //for (const ISyst * const systematic : systematicsVector)
  //systematicsToShift.push_back(systematic->ShortName());

  std::cout << "With " << sigmaShift << " sigma shift" << std::endl;

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

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};

  //std::vector<const PredictionInterp*> predictionGenerators_TRUE = {&interpGenNue_FHC_TRUE, &interpGenNue_RHC_TRUE, &interpGenNumu_FHC_TRUE, &interpGenNumu_RHC_TRUE};

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  unsigned int nThrows(N_THROWS);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  std::vector<double> chiSquaredCPCValues;
  std::vector<double> chiSquaredCPVValues;
  std::vector<double> chiSquaredValues;
  std::vector<double> deltaCPValues;

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

  for (int j = 0; j < nTestCPValues; ++j)
  {
    const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);

    //const double trueDeltaCP(0 * 3.14);

      deltaCPValues.push_back(trueDeltaCP);

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i + 1) << "/" << std::to_string(nThrows) << std::endl;

          // Get the prediction
          std::vector<Spectrum> predictionVector;
          //std::vector<Spectrum> predictionVector_M1;
          //std::vector<Spectrum> predictionVector_1;
          //std::vector<Spectrum> predictionVector_True;

	  if (MAKE_THROWS)
	  {
	    predictionVector = Get_Thrown_NO(predictionGenerators, trueDeltaCP, pot);
	  }
	  else
	  {
	    predictionVector = Get_FluxSys_NO(predictionGenerators, systematicsToShift, sigmaShift, trueDeltaCP, pot);
	    //predictionVector = Get_FluxSys_NO(predictionGenerators, systematicsToShift, 0, trueDeltaCP, pot);
	    //predictionVector_M1 = Get_FluxSys_NO(predictionGenerators, systematicsToShift, -1, trueDeltaCP, pot);
	    //predictionVector_1 = Get_FluxSys_NO(predictionGenerators, systematicsToShift, 1, trueDeltaCP, pot);
	    //predictionVector_True = Get_FluxSys_NO(predictionGenerators_TRUE, systematicsToShift, 0, trueDeltaCP, pot);
	  }

	  /*
	  predictionVector[0].ToTH1(pot)->Write("nue_0");
	  predictionVector[1].ToTH1(pot)->Write("anue_0");
	  predictionVector[2].ToTH1(pot)->Write("numu_0");
	  predictionVector[3].ToTH1(pot)->Write("anumu_0");

	  predictionVector_M1[0].ToTH1(pot)->Write("nue_M1");
	  predictionVector_M1[1].ToTH1(pot)->Write("anue_M1");
	  predictionVector_M1[2].ToTH1(pot)->Write("numu_M1");
	  predictionVector_M1[3].ToTH1(pot)->Write("anumu_M1");

	  predictionVector_1[0].ToTH1(pot)->Write("nue_1");
	  predictionVector_1[1].ToTH1(pot)->Write("anue_1");
	  predictionVector_1[2].ToTH1(pot)->Write("numu_1");
	  predictionVector_1[3].ToTH1(pot)->Write("anumu_1");

	  predictionVector_True[0].ToTH1(pot)->Write("nue_True");
	  predictionVector_True[1].ToTH1(pot)->Write("anue_True");
	  predictionVector_True[2].ToTH1(pot)->Write("numu_True");
	  predictionVector_True[3].ToTH1(pot)->Write("anumu_True");

	  return;
	  */

          // Add in reactor constraint
          // The arguments are the central value and uncertainty
          //multiExperiment.Add(new ReactorExperiment(0.088, 0.003));

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          double bestChiSquaredCPC(std::numeric_limits<float>::max());
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

	  std::map<std::string, float> bestFitPosition_CPC;
	  std::map<std::string, float> bestFitPosition_CPV;

	  std::cout << "Performing CPC fits..." << std::endl;

          // CPC Fits
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

	  //std::cout << "kFitDmSq32Scaled: " << bestFitPosition_CPC.at("kFitDmSq32Scaled") << std::endl;
	  //std::cout << "kFitSinSqTheta23: " << bestFitPosition_CPC.at("kFitSinSqTheta23") << std::endl;
	  //std::cout << "kFitSinSq2Theta12: " << bestFitPosition_CPC.at("kFitSinSq2Theta12") << std::endl;
	  //std::cout << "kFitDmSq21: " << bestFitPosition_CPC.at("kFitDmSq21") << std::endl;
	  //std::cout << "kFitTheta13: " << bestFitPosition_CPC.at("kFitTheta13") << std::endl;
	  //std::cout << "kFitRho: " << bestFitPosition_CPC.at("kFitRho") << std::endl;
	  //std::cout << "kFitDeltaInPiUnits: " << bestFitPosition_CPC.at("kFitDeltaInPiUnits") << std::endl;

	  std::cout << "Performing CPV fits..." << std::endl;

          // CPV Fits
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

	  //std::cout << "kFitDmSq32Scaled: " << bestFitPosition_CPV.at("kFitDmSq32Scaled") << std::endl;
	  //std::cout << "kFitSinSqTheta23: " << bestFitPosition_CPV.at("kFitSinSqTheta23") << std::endl;
	  //std::cout << "kFitSinSq2Theta12: " << bestFitPosition_CPV.at("kFitSinSq2Theta12") << std::endl;
	  //std::cout << "kFitDmSq21: " << bestFitPosition_CPV.at("kFitDmSq21") << std::endl;
	  //std::cout << "kFitTheta13: " << bestFitPosition_CPV.at("kFitTheta13") << std::endl;
	  //std::cout << "kFitRho: " << bestFitPosition_CPV.at("kFitRho") << std::endl;
	  //std::cout << "kFitDeltaInPiUnits: " << bestFitPosition_CPV.at("kFitDeltaInPiUnits") << std::endl;

          const float chiSquaredDifference(bestChiSquaredCPC - bestChiSquaredCPV);

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
  }

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

  tree->Fill();
  outputFile->WriteObject(tree, "tree");    
}

////////////////////////////////////////////////////////////////////////////////////////////

int GetBin(TH1D * hist, double value)
{
    for (int i = 1; i <= hist-> GetNbinsX(); ++i)
    {
        double lowEdge(hist->GetBinLowEdge(i));
        double highEdge(lowEdge + hist->GetBinWidth(i));

        if ((value > lowEdge) && (value < highEdge))
            return i;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_Thrown_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot)
{
    // Make systematic throws
    std::map<const ISyst*, double> systematicShiftMap;
    std::vector<ISyst const *> systematics(predictionGenerators.front()->GetAllSysts());

    // So, shifts have limits.. i'm pretty sure SystShifts will roll back to the boundary if we cross it
    for (const ISyst * const systematic : systematics)
        systematicShiftMap[systematic] = gRandom->Gaus();

    const SystShifts systematicShifts(systematicShiftMap);

    // Make oscillation throw
    osc::IOscCalcAdjustable* calc = DefaultOscCalc();
    //SetOscCalc_Throw_NO(calc, deltaCP);

    // Create experiments
    std::vector<Spectrum> predictionVector;

    for (const PredictionInterp *const predictionGenerator : predictionGenerators)
    {
        const Spectrum prediction(predictionGenerator->PredictSyst(calc, systematicShifts).FakeData(pot));
        predictionVector.push_back(prediction);
    }

    return predictionVector;
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_FluxSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, std::vector<std::string> &systematicsToShift, const float sigmaShift, const double deltaCP, const float pot)
{
    // Make systematic throws
    std::map<const ISyst*, double> systematicShiftMap;
    std::vector<ISyst const *> systematics(predictionGenerators.front()->GetAllSysts());

    // So, shifts have limits.. i'm pretty sure SystShifts will roll back to the boundary if we cross it
    for (const ISyst * const systematic : systematics)
    {
      std::string shortName(systematic->ShortName());

      for (std::string &sysToShiftName : systematicsToShift)
      {
          if (shortName == sysToShiftName)
          {
              std::cout << "Applying shift to shortName: " << systematic->ShortName() << std::endl;
	      systematicShiftMap[systematic] = sigmaShift;
	  }
      }
    }

    const SystShifts systematicShifts(systematicShiftMap);

    // Make oscillation calc
    osc::IOscCalcAdjustable* calc = NuFitOscCalc(1);
    calc->SetdCP(deltaCP);

    // Create experiments
    std::vector<Spectrum> predictionVector;

    for (const PredictionInterp *const predictionGenerator : predictionGenerators)
    {
        const Spectrum prediction(predictionGenerator->PredictSyst(calc, systematicShifts).AsimovData(pot));
        predictionVector.push_back(prediction);
    }

    return predictionVector;
}


////////////////////////////////////////////////////////////////////////////////////////////

Spectrum Get_CentralValue_NO(const PredictionInterp &predictionGenerator, const double deltaCP, const float pot)
{
    osc::IOscCalcAdjustable* calc = DefaultOscCalc();
    //SetOscCalc_CentralValue_NO(calc, deltaCP);

    const Spectrum prediction = predictionGenerator.Predict(calc).FakeData(pot); 

    return prediction;
}

////////////////////////////////////////////////////////////////////////////////////////////
/*
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP)
{
  calc->SetRho(2.848); 
  calc->SetDmsq21(7.39e-5);
  calc->SetDmsq32(2.451e-3);
  calc->SetTh12(0.5903);
  calc->SetTh13(0.150);
  calc->SetTh23(0.866);
  calc->SetdCP(deltaCP);
}
*/
////////////////////////////////////////////////////////////////////////////////////////////
/*
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP)
{
  //std::normal_distribution<double> th13_distribution(0.1503, 0.0023);
  std::uniform_real_distribution<double> th23_distribution(0.685, 0.886);
  std::uniform_real_distribution<double> dmsq32_distribution(2.3e-3, 2.7e-3);

  const float dmsq32 = dmsq32_distribution(generator); 
  //const float th13 = th13_distribution(generator);
  const double th23 = th23_distribution(generator);

  calc->SetRho(2.848); 
  calc->SetDmsq21(7.39e-5);
  calc->SetDmsq32(dmsq32);
  calc->SetTh12(0.5903);
  calc->SetTh13(0.150);
  calc->SetTh23(th23);
  calc->SetdCP(deltaCP);
}
*/
////////////////////////////////////////////////////////////////////////////////////////////
/*
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP)
{
  calc->SetRho(2.848); 
  calc->SetDmsq21(7.39e-5);
  calc->SetDmsq32(-2.512e-3);
  calc->SetTh12(0.5903);
  calc->SetTh13(0.151);
  calc->SetTh23(0.869);
  calc->SetdCP(deltaCP);
}
*/
////////////////////////////////////////////////////////////////////////////////////////////
/*
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP)
{
  //std::normal_distribution<double> th13_distribution(0.1510, 0.0023);
  std::uniform_real_distribution<double> th23_distribution(0.685, 0.886);
  std::uniform_real_distribution<double> dmsq32_distribution(-2.7e-3, -2.3e-3);

  const float dmsq32 = dmsq32_distribution(generator); 
  //const float th13 = th13_distribution(generator);
  const double th23 = th23_distribution(generator);

  calc->SetRho(2.848); 
  calc->SetDmsq21(7.39e-5);
  calc->SetDmsq32(dmsq32);
  calc->SetTh12(0.5903);
  calc->SetTh13(0.151);
  calc->SetTh23(th23);
  calc->SetdCP(deltaCP);
}
*/
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

  std::vector<const IFitVar*> fitVariables =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};

  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Get systematics
  std::vector<const ISyst*> systematics(predictionGenerators[0]->GetAllSysts());
  
  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});

  float chiSquared(std::numeric_limits<float>::max());

  if (FIT_SYSTEMATICS)
  {
      MinuitFitter fit(&multiExperiment, fitVariables, systematics);
      chiSquared = fit.Fit(calc)->EvalMetricVal();
  }
  else
  {
      MinuitFitter fit(&multiExperiment, fitVariables);
      chiSquared = fit.Fit(calc)->EvalMetricVal();
  }

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

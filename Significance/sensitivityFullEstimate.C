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
#include "CAFAna/Fit/StanConfig.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TRandom3.h"

#include <random>
#include <ctime>

using namespace ana;

//const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesFluxSystematicsSplitBySignCVN.root";
//const std::string INPUT_FILE_NAME = "/pnfs/dune/resilient/users/imawby/StateFiles/StateFilesAllSystematicsSplitBySign.root";
//const std::string INPUT_FILE_NAME = "CAFAna/StateFilesAllSystematicsSplitBySign.root";

void sensitivityFullEstimate(int N_THROWS, int N_TEST_DELTA_CP_VALUES, int TEST_VALUE_INDEX);
void GetSystematicsToConsider(std::vector<std::string> &systematicsToConsider);
std::vector<Spectrum> Get_Throw_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot, 
  std::map<std::string, std::vector<double>> &randomNumberMap, double &chiSquaredCPCCHEAT, double &chiSquaredCPVCHEAT, std::map<std::string, float> &bestFitPosition);
void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
  const int isHigherOctant, const int isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition, osc::IOscCalcAdjustable* bestCalcSeed, 
  SystShifts &bestSysShifts);
void PerformSeededFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, osc::IOscCalcAdjustable* calcSeed, 
  SystShifts &sysSeed, const bool isCPCFit, double &bestChiSquared, std::map<std::string, float> &bestFitPosition);
float GetBoundedGausThrowThis(float min, float max);
osc::IOscCalcAdjustable* ThrownWideOscCalcThis(int hie, std::vector<const IFitVar*> oscVars, bool flatth13, std::map<std::string, std::vector<double>> &randomNumberMap);
double GetCount(const Spectrum &spectrum);

std::default_random_engine generator;

//int N_TEST_DELTA_CP_VALUES = 100; // 200
//int N_THROWS = 500; //
const double pot = 3.5 * 1.1e21 * 40/1.13;

void sensitivityFullEstimate(int N_THROWS, int N_TEST_DELTA_CP_VALUES, int TEST_VALUE_INDEX)
{
  std::string outputFileName = "MyAnalysisOut0.root";

  std::time_t time = std::time(nullptr);
  generator.seed(time);

  // Just to register the flux systematics (HACK HACK HACK)
  DUNEFluxSystVector fluxSystematics = GetDUNEFluxSysts(30, true, false);
  
  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open("/home/epp/phrsnt/work/lblpwgtools/standardCAF/noSystematics/StateFilesNoSystematicsSplitBySign_STANDARD.root");
  //TFile * inputFile = TFile::Open("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/resilient/users/imawby/StateFiles/StateFilesAllSystematicsSplitBySign_REAL.root");
  //TFile * inputFile = TFile::Open("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/resilient/users/imawby/StateFiles/StateFilesAllSystematicsSplitBySign_STANDARD.root");
  //TFile * inputFile = TFile::Open("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/resilient/users/imawby/StateFiles/StateFilesAllSystematicsSplitBySign.root");
  //TFile * inputFile = TFile::Open("/pnfs/dune/resilient/users/imawby/StateFiles/StateFilesAllSystematicsSplitBySign.root");

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

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  TTree * tree = new TTree("tree", "tree");
  double deltaCPValues;
  std::vector<double> chiSquaredCPCValues;
  std::vector<double> chiSquaredCPVCHEATValues;
  std::vector<double> chiSquaredCPCCHEATValues;
  std::vector<double> chiSquaredCPVValues;
  std::vector<double> chiSquaredValues;

  std::vector<double> bestFitDmSq32_CPC;
  std::vector<double> bestFitTh23_CPC;
  std::vector<double> bestFitTh12_CPC;
  std::vector<double> bestFitDmSq21_CPC;
  std::vector<double> bestFitTh13_CPC;
  std::vector<double> bestFitRho_CPC;
  std::vector<double> bestFitdCP_CPC;
  std::vector<double> bestFitNueSpectrumCount_CPC;
  std::vector<double> bestFitAnueSpectrumCount_CPC;
  std::vector<double> bestFitNumuSpectrumCount_CPC;
  std::vector<double> bestFitAnumuSpectrumCount_CPC;

  std::vector<double> bestFitDmSq32_CPV;
  std::vector<double> bestFitTh23_CPV;
  std::vector<double> bestFitTh12_CPV;
  std::vector<double> bestFitDmSq21_CPV;
  std::vector<double> bestFitTh13_CPV;
  std::vector<double> bestFitRho_CPV;
  std::vector<double> bestFitdCP_CPV;
  std::vector<double> bestFitNueSpectrumCount_CPV;
  std::vector<double> bestFitAnueSpectrumCount_CPV;
  std::vector<double> bestFitNumuSpectrumCount_CPV;
  std::vector<double> bestFitAnumuSpectrumCount_CPV;

  std::vector<double> bestFitDmSq32_CHEAT;
  std::vector<double> bestFitTh23_CHEAT;
  std::vector<double> bestFitTh12_CHEAT;
  std::vector<double> bestFitDmSq21_CHEAT;
  std::vector<double> bestFitTh13_CHEAT;
  std::vector<double> bestFitRho_CHEAT;
  std::vector<double> bestFitdCP_CHEAT;

  std::map<std::string, std::vector<double>> randomNumberMap;

  tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
  tree->Branch("chiSquaredCPVValues", &chiSquaredCPVValues);
  tree->Branch("chiSquaredCPVCHEATValues", &chiSquaredCPVCHEATValues);
  tree->Branch("chiSquaredCPCCHEATValues", &chiSquaredCPCCHEATValues);
  tree->Branch("chiSquaredValues", &chiSquaredValues);
  tree->Branch("deltaCPValues", &deltaCPValues);

  tree->Branch("bestFitDmSq32_CPC", &bestFitDmSq32_CPC);
  tree->Branch("bestFitTh23_CPC", &bestFitTh23_CPC);
  tree->Branch("bestFitTh12_CPC", &bestFitTh12_CPC);
  tree->Branch("bestFitDmSq21_CPC", &bestFitDmSq21_CPC);
  tree->Branch("bestFitTh13_CPC", &bestFitTh13_CPC);
  tree->Branch("bestFitRho_CPC", &bestFitRho_CPC);
  tree->Branch("bestFitdCP_CPC", &bestFitdCP_CPC);
  tree->Branch("bestFitNueSpectrumCount_CPC", &bestFitNueSpectrumCount_CPC);
  tree->Branch("bestFitAnueSpectrumCount_CPC", &bestFitAnueSpectrumCount_CPC);
  tree->Branch("bestFitNumuSpectrumCount_CPC", &bestFitNumuSpectrumCount_CPC);
  tree->Branch("bestFitAnumuSpectrumCount_CPC", &bestFitAnumuSpectrumCount_CPC);

  tree->Branch("bestFitDmSq32_CPV", &bestFitDmSq32_CPV);
  tree->Branch("bestFitTh23_CPV", &bestFitTh23_CPV);
  tree->Branch("bestFitTh12_CPV", &bestFitTh12_CPV);
  tree->Branch("bestFitDmSq21_CPV", &bestFitDmSq21_CPV);
  tree->Branch("bestFitTh13_CPV", &bestFitTh13_CPV);
  tree->Branch("bestFitRho_CPV", &bestFitRho_CPV);
  tree->Branch("bestFitdCP_CPV", &bestFitdCP_CPV);
  tree->Branch("bestFitNueSpectrumCount_CPV", &bestFitNueSpectrumCount_CPV);
  tree->Branch("bestFitAnueSpectrumCount_CPV", &bestFitAnueSpectrumCount_CPV);
  tree->Branch("bestFitNumuSpectrumCount_CPV", &bestFitNumuSpectrumCount_CPV);
  tree->Branch("bestFitAnumuSpectrumCount_CPV", &bestFitAnumuSpectrumCount_CPV);

  tree->Branch("bestFitDmSq32_CHEAT", &bestFitDmSq32_CHEAT);
  tree->Branch("bestFitTh23_CHEAT", &bestFitTh23_CHEAT);
  tree->Branch("bestFitTh12_CHEAT", &bestFitTh12_CHEAT);
  tree->Branch("bestFitDmSq21_CHEAT", &bestFitDmSq21_CHEAT);
  tree->Branch("bestFitTh13_CHEAT", &bestFitTh13_CHEAT);
  tree->Branch("bestFitRho_CHEAT", &bestFitRho_CHEAT);
  tree->Branch("bestFitdCP_CHEAT", &bestFitdCP_CHEAT);

  // Add random numbers used in throws to the map
  std::vector<std::string> systematicsToConsider;
  GetSystematicsToConsider(systematicsToConsider);

  for (std::string &systematicToConsider : systematicsToConsider)
  {
      randomNumberMap[systematicToConsider] = std::vector<double>();
      tree->Branch(systematicToConsider.c_str(), &(randomNumberMap[systematicToConsider]));
  }

  for (std::string oscParamName : {"kFitDmSq32Scaled", "kFitSinSqTheta23", "kFitTheta13"})
  {
      randomNumberMap[oscParamName] = std::vector<double>();
      tree->Branch(oscParamName.c_str(), &(randomNumberMap[oscParamName]));
  }

  // This is left over from when i was looping over CP values...
  chiSquaredCPCValues.clear();
  chiSquaredCPVValues.clear();
  chiSquaredCPCCHEATValues.clear();
  chiSquaredCPVCHEATValues.clear();
  chiSquaredValues.clear();

  bestFitDmSq32_CPC.clear();
  bestFitTh23_CPC.clear();
  bestFitTh12_CPC.clear();
  bestFitDmSq21_CPC.clear();
  bestFitTh13_CPC.clear();
  bestFitRho_CPC.clear();
  bestFitdCP_CPC.clear();
  bestFitNueSpectrumCount_CPC.clear();
  bestFitAnueSpectrumCount_CPC.clear();
  bestFitNumuSpectrumCount_CPC.clear();
  bestFitAnumuSpectrumCount_CPC.clear();

  bestFitDmSq32_CPV.clear();
  bestFitTh23_CPV.clear();
  bestFitTh12_CPV.clear();
  bestFitDmSq21_CPV.clear();
  bestFitTh13_CPV.clear();
  bestFitRho_CPV.clear();
  bestFitdCP_CPV.clear();
  bestFitNueSpectrumCount_CPV.clear();
  bestFitAnueSpectrumCount_CPV.clear();
  bestFitNumuSpectrumCount_CPV.clear();
  bestFitAnumuSpectrumCount_CPV.clear();

  bestFitDmSq32_CHEAT.clear();
  bestFitTh23_CHEAT.clear();
  bestFitTh12_CHEAT.clear();
  bestFitDmSq21_CHEAT.clear();
  bestFitTh13_CHEAT.clear();
  bestFitRho_CHEAT.clear();
  bestFitdCP_CHEAT.clear();

  for (auto &entry : randomNumberMap)
      entry.second.clear();

  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  unsigned int nThrows(N_THROWS);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  const double trueDeltaCP(static_cast<float>(TEST_VALUE_INDEX) * stepSizeCP);

  deltaCPValues = trueDeltaCP;

  //Perform throws
      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "nThrow: " << std::to_string(i + 1) << "/" << std::to_string(nThrows) << std::endl;

          double bestChiSquareCPCCHEAT(std::numeric_limits<float>::max());
          double bestChiSquareCPVCHEAT(std::numeric_limits<float>::max());
          std::map<std::string, float> bestFitPosition_CHEAT;

          // Get the prediction (we need a cheatCalcSeed and a cheatSysSeed)
          std::vector<Spectrum> predictionVector;
	  predictionVector = Get_Throw_NO(predictionGenerators, trueDeltaCP, pot, randomNumberMap, bestChiSquareCPCCHEAT, bestChiSquareCPVCHEAT, bestFitPosition_CHEAT);

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          double bestChiSquaredCPC(std::numeric_limits<float>::max());
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

          std::map<std::string, float> bestFitPosition_CPC;
          std::map<std::string, float> bestFitPosition_CPV;

          osc::IOscCalcAdjustable* bestCalcSeed = NuFitOscCalc(1);
          SystShifts bestSysSeed;
      
          // CPC Fits
          std::cout << "Performing CPC fits..." << std::endl;
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), 1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), 1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), -1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), -1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
      
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), 1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), 1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), -1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), -1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);

          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), 1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), 1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), -1, 1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), -1, -1, true, bestChiSquaredCPC, bestFitPosition_CPC, bestCalcSeed, bestSysSeed);

          // CPV Fits (These guys cannot touch the best seed objects - sorry that i'm pure trash)
          std::cout << "Performing CPV fits..." << std::endl;
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), 1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), 1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), -1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.0 * TMath::Pi(), -1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);

          PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), 1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), 1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), -1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 0.5 * TMath::Pi(), -1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);

          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), 1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), 1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), -1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.0 * TMath::Pi(), -1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);

          PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), 1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), 1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), -1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 1.5 * TMath::Pi(), -1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);

          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), 1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), 1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), -1, 1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
          PerformFit(predictionGenerators, predictionVector, 2.0 * TMath::Pi(), -1, -1, false, bestChiSquaredCPV, bestFitPosition_CPV, bestCalcSeed, bestSysSeed);
      
          // Okay, now to make sure that the CPV fit finds a better position compared to the CPC fit (might have to seed it at cheated position) UGGGHH THESE DONT STORE THE BEST FIT POSITION. -.-
          std::cout << "Performing CPC Seeded Fit..." << std::endl;
          PerformSeededFit(predictionGenerators, predictionVector, bestCalcSeed, bestSysSeed, false, bestChiSquaredCPV, bestFitPosition_CPV);

          // Use the cheated seed as well to solve the asymmetry problem
          const float chiSquaredDifference(std::min(bestChiSquaredCPC, bestChiSquareCPCCHEAT) - std::min(bestChiSquaredCPV, bestChiSquareCPVCHEAT));

          //std::cout << "////////////////////////////" << std::endl;
          //std::cout << "trueDeltaCP: " << (trueDeltaCP * 180 / 3.14) << std::endl;
          //std::cout << "bestChiSquaredCPC: " << bestChiSquaredCPC << std::endl;
          //std::cout << "bestChiSquaredCPV: " << bestChiSquaredCPV << std::endl;
          //std::cout << "chiSquaredDifference: " << chiSquaredDifference << std::endl;
          //std::cout << "significance: " << TMath::Sqrt(std::fabs(chiSquaredDifference)) << std::endl;
          //std::cout << "////////////////////////////" << std::endl;

          chiSquaredCPCValues.push_back(bestChiSquaredCPC);
          chiSquaredCPVValues.push_back(bestChiSquaredCPV);
          chiSquaredCPCCHEATValues.push_back(bestChiSquareCPCCHEAT);
          chiSquaredCPVCHEATValues.push_back(bestChiSquareCPVCHEAT);
          chiSquaredValues.push_back(chiSquaredDifference);

          /*
          bestFitDmSq32_CPC.push_back(bestFitPosition_CPC.find("kFitDmSq32Scaled") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq32Scaled") : -999.0);
          bestFitTh23_CPC.push_back(bestFitPosition_CPC.find("kFitSinSqTheta23") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSqTheta23") : -999.0);
          bestFitTh12_CPC.push_back(bestFitPosition_CPC.find("kFitSinSq2Theta12") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitSinSq2Theta12") : -999.0);
          bestFitDmSq21_CPC.push_back(bestFitPosition_CPC.find("kFitDmSq21") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDmSq21") : -999.0);
          bestFitTh13_CPC.push_back(bestFitPosition_CPC.find("kFitTheta13") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitTheta13") : -999.0);
          bestFitRho_CPC.push_back(bestFitPosition_CPC.find("kFitRho") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitRho") : -999.0);
          bestFitdCP_CPC.push_back(bestFitPosition_CPC.find("kFitDeltaInPiUnits") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("kFitDeltaInPiUnits") : -999.0);
          */
          bestFitDmSq32_CPC.push_back(-999.0);
          bestFitTh23_CPC.push_back(-999.0);
          bestFitTh12_CPC.push_back(-999.0);
          bestFitDmSq21_CPC.push_back(-999.0);
          bestFitTh13_CPC.push_back(-999.0);
          bestFitRho_CPC.push_back(-999.0);
          bestFitdCP_CPC.push_back(-999.0);

          bestFitNueSpectrumCount_CPC.push_back(bestFitPosition_CPC.find("NueSpectrumCount") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("NueSpectrumCount") : -999.0);
          bestFitAnueSpectrumCount_CPC.push_back(bestFitPosition_CPC.find("AnueSpectrumCount") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("AnueSpectrumCount") : -999.0);
          bestFitNumuSpectrumCount_CPC.push_back(bestFitPosition_CPC.find("NumuSpectrumCount") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("NumuSpectrumCount") : -999.0);
          bestFitAnumuSpectrumCount_CPC.push_back(bestFitPosition_CPC.find("AnumuSpectrumCount") != bestFitPosition_CPC.end() ? bestFitPosition_CPC.at("AnumuSpectrumCount") : -999.0);

          bestFitDmSq32_CPV.push_back(bestFitPosition_CPV.find("kFitDmSq32Scaled") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDmSq32Scaled") : -999.0);
          bestFitTh23_CPV.push_back(bestFitPosition_CPV.find("kFitSinSqTheta23") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitSinSqTheta23") : -999.0);
          bestFitTh12_CPV.push_back(bestFitPosition_CPV.find("kFitSinSq2Theta12") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitSinSq2Theta12") : -999.0);
          bestFitDmSq21_CPV.push_back(bestFitPosition_CPV.find("kFitDmSq21") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDmSq21") : -999.0);
          bestFitTh13_CPV.push_back(bestFitPosition_CPV.find("kFitTheta13") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitTheta13") : -999.0);
          bestFitRho_CPV.push_back(bestFitPosition_CPV.find("kFitRho") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitRho") : -999.0);
          bestFitdCP_CPV.push_back(bestFitPosition_CPV.find("kFitDeltaInPiUnits") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("kFitDeltaInPiUnits") : -999.0);
          bestFitNueSpectrumCount_CPV.push_back(bestFitPosition_CPV.find("NueSpectrumCount") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("NueSpectrumCount") : -999.0);
          bestFitAnueSpectrumCount_CPV.push_back(bestFitPosition_CPV.find("AnueSpectrumCount") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("AnueSpectrumCount") : -999.0);
          bestFitNumuSpectrumCount_CPV.push_back(bestFitPosition_CPV.find("NumuSpectrumCount") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("NumuSpectrumCount") : -999.0);
          bestFitAnumuSpectrumCount_CPV.push_back(bestFitPosition_CPV.find("AnumuSpectrumCount") != bestFitPosition_CPV.end() ? bestFitPosition_CPV.at("AnumuSpectrumCount") : -999.0);

          bestFitDmSq32_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitDmSq32Scaled") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitDmSq32Scaled") : -999.0);
          bestFitTh23_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitSinSqTheta23") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitSinSqTheta23") : -999.0);
          bestFitTh12_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitSinSq2Theta12") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitSinSq2Theta12") : -999.0);
          bestFitDmSq21_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitDmSq21") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitDmSq21") : -999.0);
          bestFitTh13_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitTheta13") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitTheta13") : -999.0);
          bestFitRho_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitRho") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitRho") : -999.0);
          bestFitdCP_CHEAT.push_back(bestFitPosition_CHEAT.find("kFitDeltaInPiUnits") != bestFitPosition_CHEAT.end() ? bestFitPosition_CHEAT.at("kFitDeltaInPiUnits") : -999.0);
      }  

      tree->Fill();

  outputFile->WriteObject(tree, "tree");    
  outputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////

void GetSystematicsToConsider(std::vector<std::string> &systematicsToConsider)
{
  systematicsToConsider = {
    "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt", 
    "ChargedHadUncorrFD", "UncorrFDHadSqrt", "UncorrFDHadInvSqrt", "ChargedHadResFD",
    "EScaleMuLArFD", "UncorrFDMuSqrt", "UncorrFDMuInvSqrt", "MuonResFD",
    "NUncorrFD", "UncorrFDNSqrt", "UncorrFDNInvSqrt", "NResFD",
    "EMUncorrFD", "UncorrFDEMSqrt", "UncorrFDEMInvSqrt", "EMResFD"};

  DUNEFluxSystVector fluxSystematics = GetDUNEFluxSysts(30, true, false);

  for (const ISyst * const fluxSystematic : fluxSystematics)
      systematicsToConsider.push_back(fluxSystematic->ShortName());

  std::vector<ana::ISyst const *> xsecSystematics = GetListOfSysts(false, true, false, false, true, false, true, false, 0, false);

  for (const ISyst * const xsecSystematic : xsecSystematics)
      systematicsToConsider.push_back(xsecSystematic->ShortName());
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_Throw_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot, 
  std::map<std::string, std::vector<double>> &randomNumberMap, double &chiSquaredCPCCheat, double &chiSquaredCPVCheat, std::map<std::string, float> &bestFitPosition)
{
    std::vector<std::string> systematicsToConsider;
    GetSystematicsToConsider(systematicsToConsider);

    // Make systematic throws
    SystShifts systematicShifts;

    /*
    std::vector<ISyst const *> systematics(predictionGenerators.front()->GetAllSysts());

    for (const ISyst * const systematic : systematics)
    {
        std::string shortName(systematic->ShortName());

        for (std::string &sysToShiftName : systematicsToConsider)
        {
            if (shortName == sysToShiftName)
            {
                const float shift = GetBoundedGausThrowThis(systematic->Min() * 0.8, systematic->Max() * 0.8);
                randomNumberMap[shortName].push_back(shift);
                systematicShifts.SetShift(systematic, shift);
            }
        }
    }
    */
    // Make oscillation calc throw (1 is for NO)
    std::vector<const IFitVar *> oscVars = GetOscVars("dmsq32:th23:th13", 1);
    osc::IOscCalcAdjustable* calc = ThrownWideOscCalcThis(1, oscVars, false, randomNumberMap);
    calc->SetdCP(deltaCP);

    // Create experiments (add in statistical fluctuations)
    std::vector<Spectrum> predictionVector;
    for (const PredictionInterp *const predictionGenerator : predictionGenerators)
    {
      //const Spectrum prediction(predictionGenerator->PredictSyst(calc, systematicShifts).MockData(pot));
        const Spectrum prediction(predictionGenerator->Predict(calc).AsimovData(pot));
        predictionVector.push_back(prediction);
    }

    /*
  std::cout << "----- calc seed before -------" << std::endl;

  std::cout << "rho: " << calc->GetRho() << std::endl;
  std::cout << "dmsq21: " << calc->GetDmsq21() << std::endl;
  std::cout << "dmsq32: " << calc->GetDmsq32() << std::endl;
  std::cout << "theta12: " << calc->GetTh12() << std::endl;
  std::cout << "theta13: " << calc->GetTh13() << std::endl;
  std::cout << "theta23: " << calc->GetTh23() << std::endl;
  std::cout << "DeltaCP: " << calc->GetdCP() << std::endl;
    */

    // For the seeded CPC fit i need to make copies
    // Because the calc and syshift object will change -.- 
    SystShifts systematicShiftsZero;
    SystShifts systematicShiftsPi;
    /*
    auto activeSys = systematicShifts.ActiveSysts();
    for (auto &entry : activeSys)
    {
        systematicShiftsZero.SetShift(entry, systematicShifts.GetShift(entry));
        systematicShiftsPi.SetShift(entry, systematicShifts.GetShift(entry));
    }
    */
    osc::IOscCalcAdjustable* calcPi = NuFitOscCalc(1);
    calcPi->SetRho(calc->GetRho());
    calcPi->SetDmsq21(calc->GetDmsq21());
    calcPi->SetDmsq32(calc->GetDmsq32());
    calcPi->SetTh12(calc->GetTh12());
    calcPi->SetTh13(calc->GetTh13());
    calcPi->SetTh23(calc->GetTh23());
    calcPi->SetdCP(1.0 * TMath::Pi());

    osc::IOscCalcAdjustable* calcZero = NuFitOscCalc(1);
    calcZero->SetRho(calc->GetRho());
    calcZero->SetDmsq21(calc->GetDmsq21());
    calcZero->SetDmsq32(calc->GetDmsq32());
    calcZero->SetTh12(calc->GetTh12());
    calcZero->SetTh13(calc->GetTh13());
    calcZero->SetTh23(calc->GetTh23());
    calcZero->SetdCP(0.0 * TMath::Pi());
 
    // Perform fits to check the true position
    std::cout << "Performing Truth CPC Seeded Fit..." << std::endl;
    std::map<std::string, float> dummyBestFitPosition;
    double chiSquaredZero(std::numeric_limits<float>::max()), chiSquaredPi(std::numeric_limits<float>::max());
    PerformSeededFit(predictionGenerators, predictionVector, calcZero, systematicShiftsZero, true, chiSquaredZero, dummyBestFitPosition);
    PerformSeededFit(predictionGenerators, predictionVector, calcPi, systematicShiftsPi, true, chiSquaredPi, dummyBestFitPosition);
    chiSquaredCPCCheat = std::min(chiSquaredZero, chiSquaredPi);

    // Set the CPC cheated seed
    osc::IOscCalcAdjustable* cheatCalcSeedCPC = NuFitOscCalc(1);
    SystShifts cheatSysSeedCPC;

    if (chiSquaredZero < chiSquaredPi)
    {
        cheatCalcSeedCPC->SetRho(calcZero->GetRho());
        cheatCalcSeedCPC->SetDmsq21(calcZero->GetDmsq21());
        cheatCalcSeedCPC->SetDmsq32(calcZero->GetDmsq32());
        cheatCalcSeedCPC->SetTh12(calcZero->GetTh12());
        cheatCalcSeedCPC->SetTh13(calcZero->GetTh13());
        cheatCalcSeedCPC->SetTh23(calcZero->GetTh23());
        cheatCalcSeedCPC->SetdCP(calcZero->GetdCP());

	/*
        for (auto &entry : activeSys)
          cheatSysSeedCPC.SetShift(entry, systematicShiftsZero.GetShift(entry));
	*/
    }
    else
    {
        cheatCalcSeedCPC->SetRho(calcPi->GetRho());
        cheatCalcSeedCPC->SetDmsq21(calcPi->GetDmsq21());
        cheatCalcSeedCPC->SetDmsq32(calcPi->GetDmsq32());
        cheatCalcSeedCPC->SetTh12(calcPi->GetTh12());
        cheatCalcSeedCPC->SetTh13(calcPi->GetTh13());
        cheatCalcSeedCPC->SetTh23(calcPi->GetTh23());
        cheatCalcSeedCPC->SetdCP(calcPi->GetdCP());

	/*
        for (auto &entry : activeSys)
          cheatSysSeedCPC.SetShift(entry, systematicShiftsPi.GetShift(entry));
	*/
    }

    std::cout << "Performing Truth CPV Seeded Fit..." << std::endl;
    PerformSeededFit(predictionGenerators, predictionVector, calc, systematicShifts, false, chiSquaredCPVCheat, bestFitPosition);

    // Also seed at the CPC position, to prevent negatives.. lol i hate this.
    PerformSeededFit(predictionGenerators, predictionVector, cheatCalcSeedCPC, cheatSysSeedCPC, false, chiSquaredCPVCheat, bestFitPosition);

    return predictionVector;
}

////////////////////////////////////////////////////////////////////////////////////////////

float GetBoundedGausThrowThis(float min, float max)
{
   std::normal_distribution<double> gausDistribution(0, 1);

   float shift = -9999.0;

   while ((shift < min) || (shift > max))
     shift = gausDistribution(generator);

   return shift;
}

////////////////////////////////////////////////////////////////////////////////////////////

osc::IOscCalcAdjustable* ThrownWideOscCalcThis(int hie, std::vector<const IFitVar*> oscVars, bool flatth13, std::map<std::string, std::vector<double>> &randomNumberMap)
{
  std::uniform_real_distribution<double> uniformDistributionDmSq23(2.3e-3, 2.7e-3);
  std::uniform_real_distribution<double> uniformDistributionTh23(asin(sqrt(0.4)), asin(sqrt(0.6)));
  std::normal_distribution<double> gausDistribution(0, 1);

  const double kNuFitTh13CVNH = 8.61 * TMath::Pi()/180;
  const double kNuFitTh13CVIH = 8.65 * TMath::Pi()/180;
  const double kNuFitTh13ErrNH = ((8.99-8.22)/6) * TMath::Pi()/180;
  const double kNuFitTh13ErrIH = ((9.03-8.27)/6) * TMath::Pi()/180;

    assert(hie == +1 || hie == -1);

    osc::IOscCalcAdjustable* ret = NuFitOscCalc(hie);

    // Throw 12 and rho within errors if here
    //if (HasVar(oscVars, kFitRho.ShortName()))
    //ret->SetRho(kEarthDensity*(1+0.02*gRandom->Gaus()));

    //if (HasVar(oscVars, kFitDmSq21.ShortName()) or
	//HasVar(oscVars, kFitDmSq21Scaled.ShortName()))
    //ret->SetDmsq21(kNuFitDmsq21CV+kNuFitDmsq21Err*gRandom->Gaus());

    //if (HasVar(oscVars, kFitSinSq2Theta12.ShortName()))
    //ret->SetTh12(kNuFitTh12CV+kNuFitTh12Err*gRandom->Gaus());

    // Throw dmsq32 flat between 2.3 and 2.7 in the correct hierarchy
    if (HasVar(oscVars, kFitDmSq32Scaled.ShortName()) or
	HasVar(oscVars, kFitDmSq32NHScaled.ShortName()) or
	HasVar(oscVars, kFitDmSq32IHScaled.ShortName()))
    {
        const double shift(uniformDistributionDmSq23(generator));
        ret->SetDmsq32(float(hie)*shift);
        randomNumberMap["kFitDmSq32Scaled"].push_back(shift);
    }

    // Throw sin2th23 flat between 0.4 and 0.6
    if (HasVar(oscVars, kFitSinSqTheta23.ShortName()))
    {
        const double shift(uniformDistributionTh23(generator));
        ret->SetTh23(shift);
        randomNumberMap["kFitSinSqTheta23"].push_back(shift);
    }
    
    //if (HasVar(oscVars, kFitSinSqTheta23UpperOctant.ShortName()))
    //ret->SetTh23(gRandom->Uniform(th23_med, th23_high));

    //if (HasVar(oscVars, kFitSinSqTheta23LowerOctant.ShortName()))
    //ret->SetTh23(gRandom->Uniform(th23_low, th23_med));

    // Throw th13
    if (HasVar(oscVars, kFitTheta13.ShortName()))
    {
        //if (flatth13)
        //ret->SetTh13(gRandom->Uniform(0.13, 0.2));
        //else
        const double shift(gausDistribution(generator));
	    ret->SetTh13((hie > 0 ? kNuFitTh13CVNH : kNuFitTh13CVIH) *
          (1 + (hie > 0 ? kNuFitTh13ErrNH : kNuFitTh13ErrIH)*shift));
        randomNumberMap["kFitTheta13"].push_back(shift);
    }

    // dCP is just flat
    //if (HasVar(oscVars, kFitDeltaInPiUnits.ShortName()))
    //ret->SetdCP(gRandom->Uniform(-1*TMath::Pi(), TMath::Pi()));

    return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////

void PerformFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, const double deltaCPSeed, 
  const int isHigherOctant, const int isPositiveHierarchy, bool fitCPC, double &bestChiSquared, std::map<std::string, float> &bestFitPosition, osc::IOscCalcAdjustable* bestCalcSeed, 
  SystShifts &bestSysShifts)
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
  /*
  std::vector<const ISyst*> systematicsToFit;
  std::vector<std::string> systematicsToConsider;
  GetSystematicsToConsider(systematicsToConsider);

  std::vector<ISyst const *> systematics(predictionGenerators.front()->GetAllSysts());

  for (const ISyst * const systematic : systematics)
  {
      std::string shortName(systematic->ShortName());

      for (std::string &sysToShiftName : systematicsToConsider)
      {
          if (shortName == sysToShiftName)
              systematicsToFit.push_back(systematic);
      }
  }
  */
  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  // Create multi-experiment
  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});

  // Perform fit, retain best fit position
  SystShifts fitSys;
  float chiSquared(std::numeric_limits<float>::max());

  MinuitFitter fit(&multiExperiment, fitVariables);//, systematicsToFit);
  //const std::unique_ptr<IFitter::IFitSummary> fitSummary(fit.Fit(calc, fitSys, SeedList(), {}, ana::IFitter::kQuiet));
  //chiSquared = fitSummary->EvalMetricVal();

  chiSquared = fit.Fit(calc)->EvalMetricVal();

  //bool isGood(fitSummary->IsFitGood()); //idk what to do with this..
 
  if (chiSquared < bestChiSquared)
  {
       bestChiSquared = chiSquared;

       bestFitPosition["kFitDmSq32Scaled"] = calc->GetDmsq32();
       bestFitPosition["kFitSinSqTheta23"] = calc->GetTh23();
       bestFitPosition["kFitSinSq2Theta12"] = calc->GetTh12();
       bestFitPosition["kFitDmSq21"] = calc->GetDmsq21();
       bestFitPosition["kFitTheta13"] = calc->GetTh13();
       bestFitPosition["kFitRho"] = calc->GetRho();
       bestFitPosition["kFitDeltaInPiUnits"] = calc->GetdCP();

       if (fitCPC)  
       {
	 /*
           // Set the systematic seed...
           bestSysShifts.ResetToNominal();
           std::vector<const ISyst*> activeSys = fitSys.ActiveSysts();

           for (auto entry : activeSys)
               bestSysShifts.SetShift(entry, fitSys.GetShift(entry));
	 */

           // Set the calculator seed
           bestCalcSeed->SetDmsq32(calc->GetDmsq32());
           bestCalcSeed->SetTh23(calc->GetTh23());
           bestCalcSeed->SetTh12(calc->GetTh12());
           bestCalcSeed->SetDmsq21(calc->GetDmsq21());
           bestCalcSeed->SetTh13(calc->GetTh13());
           bestCalcSeed->SetRho(calc->GetRho());
           bestCalcSeed->SetdCP(calc->GetdCP());
       }

       //const Spectrum expectationNue(predictionGenerators[0]->PredictSyst(calc, fitSys).AsimovData(pot));
       //const Spectrum expectationAnue(predictionGenerators[1]->PredictSyst(calc, fitSys).AsimovData(pot));
       //const Spectrum expectationNumu(predictionGenerators[2]->PredictSyst(calc, fitSys).AsimovData(pot));
       //const Spectrum expectationAnumu(predictionGenerators[3]->PredictSyst(calc, fitSys).AsimovData(pot));

       const Spectrum expectationNue(predictionGenerators[0]->Predict(calc).AsimovData(pot));
       const Spectrum expectationAnue(predictionGenerators[1]->Predict(calc).AsimovData(pot));
       const Spectrum expectationNumu(predictionGenerators[2]->Predict(calc).AsimovData(pot));
       const Spectrum expectationAnumu(predictionGenerators[3]->Predict(calc).AsimovData(pot));

       bestFitPosition["NueSpectrumCount"] = GetCount(expectationNue);
       bestFitPosition["AnueSpectrumCount"] = GetCount(expectationAnue);
       bestFitPosition["NumuSpectrumCount"] = GetCount(expectationNumu);
       bestFitPosition["AnumuSpectrumCount"] = GetCount(expectationAnumu);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////

void PerformSeededFit(std::vector<const PredictionInterp*> &predictionGenerators, const std::vector<Spectrum> &predictionVector, osc::IOscCalcAdjustable* calcSeed, 
  SystShifts &sysSeed, const bool isCPCFit, double &bestChiSquared, std::map<std::string, float> &bestFitPosition)
{
  // Create the oscillation seed
  std::map<const IFitVar*, double> oscSeedMap;
  oscSeedMap[&kFitDmSq32Scaled] = kFitDmSq32Scaled.GetValue(calcSeed);
  oscSeedMap[&kFitSinSqTheta23] = kFitSinSqTheta23.GetValue(calcSeed);
  oscSeedMap[&kFitSinSq2Theta12] = kFitSinSq2Theta12.GetValue(calcSeed);
  oscSeedMap[&kFitDmSq21] = kFitDmSq21.GetValue(calcSeed);
  oscSeedMap[&kFitTheta13] = kFitTheta13.GetValue(calcSeed);
  oscSeedMap[&kFitRho] = kFitRho.GetValue(calcSeed);

  if (!isCPCFit)
      oscSeedMap[&kFitDeltaInPiUnits] = kFitDeltaInPiUnits.GetValue(calcSeed);

  const Seed oscSeed(oscSeedMap);
  const SeedList oscSeedList({oscSeed});

  // Get fit variables
  std::vector<const IFitVar*> fitVariables =
    {&kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitSinSq2Theta12, &kFitDmSq21,
     &kFitTheta13, &kFitRho};

  if (!isCPCFit)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  // Get systematics
  std::vector<const ISyst*> systematicsToFit;
  std::vector<std::string> systematicsToConsider;
  /*
  GetSystematicsToConsider(systematicsToConsider);

  std::vector<ISyst const *> systematics(predictionGenerators.front()->GetAllSysts());

  for (const ISyst * const systematic : systematics)
  {
      std::string shortName(systematic->ShortName());

      for (std::string &sysToShiftName : systematicsToConsider)
      {
          if (shortName == sysToShiftName)
	      systematicsToFit.push_back(systematic);
      }
  }
  */

  // Turn into an experiment (I agree, it is annoying they're here.. but be here they must because the penalty is an experiment)
  const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
  const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
  const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
  const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);

  // include the th13, rho, th12, dm21 penalty.. 
  int hie = calcSeed->GetDmsq32() > 0 ? 1 : -1;
  int oct = std::pow(sin(calcSeed->GetTh23()), 2.0) > 0.5 ? 1 : -1;

  Penalizer_GlbLike penalty(hie, oct, true, false, false, 0);

  // Define multi-experiment
  MultiExperiment multiExperiment({&experimentNueFHC, &experimentNueRHC, &experimentNumuFHC, &experimentNumuRHC, &penalty});

  // Perform fit, fitCalc and fitSys will catch the fit position
  SystShifts fitSys;

  MinuitFitter fit(&multiExperiment, fitVariables);//, systematicsToFit);
  //float chiSquared = fit.Fit(calcSeed, fitSys, oscSeedList, {sysSeed})->EvalMetricVal();
  float chiSquared = fit.Fit(calcSeed)->EvalMetricVal();

  if (chiSquared < bestChiSquared)
  {
      bestChiSquared = chiSquared;

      bestFitPosition["kFitDmSq32Scaled"] = calcSeed->GetDmsq32();
      bestFitPosition["kFitSinSqTheta23"] = calcSeed->GetTh23();
      bestFitPosition["kFitSinSq2Theta12"] = calcSeed->GetTh12();
      bestFitPosition["kFitDmSq21"] = calcSeed->GetDmsq21();
      bestFitPosition["kFitTheta13"] = calcSeed->GetTh13();
      bestFitPosition["kFitRho"] = calcSeed->GetRho();
      bestFitPosition["kFitDeltaInPiUnits"] = calcSeed->GetdCP();
  }

  // could change them at the end?? (the calc will change...)
  /*
  sysSeed.ResetToNominal();
  std::vector<const ISyst*> activeSys = fitSys.ActiveSysts();

  for (auto entry : activeSys)
      sysSeed.SetShift(entry, fitSys.GetShift(entry));
  */
}

////////////////////////////////////////////////////////////////////////////////////////////

double GetCount(const Spectrum &spectrum)
{
    double totalCount = 0;

    TH1 * hist = spectrum.ToTH1(pot);
    int nBins = hist->GetNbinsX();

    for (int i = 1; i <=nBins; ++i)
    {
        const double binContent = hist->GetBinContent(i);
        totalCount += binContent;
    }

    return totalCount;
}

  /*std::cout << "after: " << std::endl;
  std::cout << "rho: " << calc->GetRho() << std::endl;
  std::cout << "dmsq21: " << calc->GetDmsq21() << std::endl;
  std::cout << "dmsq32: " << calc->GetDmsq32() << std::endl;
  std::cout << "theta12: " << calc->GetTh12() << std::endl;
  std::cout << "theta13: " << calc->GetTh13() << std::endl;
  std::cout << "theta23: " << calc->GetTh23() << std::endl;
  std::cout << "DeltaCP: " << calc->GetdCP() << std::endl;*/

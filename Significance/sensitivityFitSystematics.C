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
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Prediction/PredictionInterp.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TRandom3.h"

#include <random>

using namespace ana;

const std::string INPUT_FILE_NAME = "/dune/data/users/imawby/standardCAF/StateFilesDetectorSystematics.root";

void PerformFit(MultiExperiment &experiment, std::vector<const ISyst*> &systematics, const double deltaCPSeed, const bool isLowerOctant, bool fitCPC, double &bestChiSquared);
std::vector<Spectrum> Get_Thrown_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);
Spectrum Get_CentralValue_NO(const PredictionInterp &predictionGenerator, const double deltaCP, const float pot);
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

std::default_random_engine generator;

int N_TEST_DELTA_CP_VALUES = 2; // 10
bool MAKE_THROWS = false;
int N_THROWS = MAKE_THROWS ? 200 : 1; // 500
bool FIT_SYSTEMATICS = false;

void sensitivityFitSystematics()
{
  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  PredictionInterp& interpGenNue_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_IZZLE").release();
  //PredictionInterp& interpGenNue_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_IZZLE").release();
  //PredictionInterp& interpGenNumu_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_IZZLE").release();
  //PredictionInterp& interpGenNumu_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_IZZLE").release();

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE}; //, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile("/dune/data/users/imawby/standardCAF/IZZLESensitivityPlotsSystematics_NO.root", "CREATE");
  std::cout << "created output file" << std::endl;


  // Evaluate each test value
  int nTestCPValues(N_TEST_DELTA_CP_VALUES);
  unsigned int nThrows(N_THROWS);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  TTree * tree = new TTree("tree", "tree");
  std::vector<std::vector<double>> chiSquaredValues(nTestCPValues);
  std::vector<double> deltaCPValues(nTestCPValues);

  for (int j = 0; j < nTestCPValues; ++j)
  {
      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);
      deltaCPValues[j] = trueDeltaCP;

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i) << "/" << std::to_string(nThrows) << std::endl;

          // Get the prediction
          std::vector<Spectrum> predictionVector(Get_Thrown_NO(predictionGenerators, trueDeltaCP, pot));

          // Turn into an experiment
          MultiExperiment multiExperiment;
          //for (unsigned int i = 0; i < predictionGenerators.size(); ++i)
          // {
              const SingleSampleExperiment experiment(predictionGenerators[1], predictionVector[1]);
              multiExperiment.Add(&experiment);
              //}

          // Add in reactor constraint
          // The arguments are the central value and uncertainty
          multiExperiment.Add(new ReactorExperiment(0.088, 0.003));

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          std::vector<const ISyst*> systematics(interpGenNue_FHC_IZZLE.GetAllSysts());
          double bestChiSquaredCPC(std::numeric_limits<float>::max());
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

          // CPC Fits
          PerformFit(multiExperiment, systematics, 0.0 * TMath::Pi(), true, true, bestChiSquaredCPC);
          PerformFit(multiExperiment, systematics, 0.0 * TMath::Pi(), false, true, bestChiSquaredCPC);

          PerformFit(multiExperiment, systematics, 1.0 * TMath::Pi(), true, true, bestChiSquaredCPC);
          PerformFit(multiExperiment, systematics, 1.0 * TMath::Pi(), false, true, bestChiSquaredCPC);

          PerformFit(multiExperiment, systematics, 2.0 * TMath::Pi(), true, true, bestChiSquaredCPC);
          PerformFit(multiExperiment, systematics, 2.0 * TMath::Pi(), false, true, bestChiSquaredCPC);

          // CPV Fits
          PerformFit(multiExperiment, systematics, 0.0 * TMath::Pi(), true, false, bestChiSquaredCPV);
          PerformFit(multiExperiment, systematics, 0.0 * TMath::Pi(), false, false, bestChiSquaredCPV);

          PerformFit(multiExperiment, systematics, 1.0 * TMath::Pi(), true, false, bestChiSquaredCPV);
          PerformFit(multiExperiment, systematics, 1.0 * TMath::Pi(), false, false, bestChiSquaredCPV);

          PerformFit(multiExperiment, systematics, 2.0 * TMath::Pi(), true, false, bestChiSquaredCPV);
          PerformFit(multiExperiment, systematics, 2.0 * TMath::Pi(), false, false, bestChiSquaredCPV);

          const float chiSquaredDifference(bestChiSquaredCPC - bestChiSquaredCPV);

          std::cout << "////////////////////////////" << std::endl;
          std::cout << "trueDeltaCP: " << (trueDeltaCP * 180 / 3.14) << std::endl;
          std::cout << "bestChiSquaredCPC: " << bestChiSquaredCPC << std::endl;
          std::cout << "bestChiSquaredCPV: " << bestChiSquaredCPV << std::endl;
          std::cout << "chiSquaredDifference: " << chiSquaredDifference << std::endl;
          std::cout << "significance: " << TMath::Sqrt(chiSquaredDifference) << std::endl;
          std::cout << "////////////////////////////" << std::endl;

          chiSquaredValues[j].push_back(bestChiSquaredCPC);
      }

      tree->Branch(("chiSquaredValues_" + std::to_string(j)).c_str(), &chiSquaredValues[j]);
      tree->Branch(("deltaCPValues_" + std::to_string(j)).c_str(), &deltaCPValues[j]);  
  }

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
    SetOscCalc_Throw_NO(calc, deltaCP);

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

Spectrum Get_CentralValue_NO(const PredictionInterp &predictionGenerator, const double deltaCP, const float pot)
{
    osc::IOscCalcAdjustable* calc = DefaultOscCalc();
    SetOscCalc_CentralValue_NO(calc, deltaCP);

    const Spectrum prediction = predictionGenerator.Predict(calc).FakeData(pot); 

    return prediction;
}

////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////

void PerformFit(MultiExperiment &experiment, std::vector<const ISyst*> &systematics, const double deltaCPSeed, const bool isLowerOctant, bool fitCPC, double &bestChiSquared)
{
  // Define the fit variables
    std::vector<const IFitVar *> fitVariables = {&kFitDmSq32NHScaled, &kFitSinSq2Theta12, &kFitDmSq21, &kFitRho, &kFitTheta13};

  isLowerOctant ? fitVariables.push_back(&kFitSinSqTheta23LowerOctant) : fitVariables.push_back(&kFitSinSqTheta23UpperOctant);

  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  /* this is a longer list of things i should be fitting
  std::vector<const IFitVar*> oscVars =
      {&kFitDmSq32Scaled, &kFitSinSqTheta23,
       &kFitSinSq2Theta12, &kFitDmSq21,
       &kFitTheta13, &kFitRho};
  */

  // Set the oscillation calculator seeed
  osc::IOscCalcAdjustable* calc = DefaultOscCalc();
  SetOscCalc_CentralValue_NO(calc, deltaCPSeed);

  if (isLowerOctant)
      calc->SetTh23(0.735);
  else
      calc->SetTh23(0.835);


  float chiSquared(std::numeric_limits<float>::max());

  if (FIT_SYSTEMATICS)
  {
      MinuitFitter fit(&experiment, fitVariables, systematics);
      chiSquared = fit.Fit(calc)->EvalMetricVal();
  }
  else
  {
      MinuitFitter fit(&experiment, fitVariables);
      chiSquared = fit.Fit(calc)->EvalMetricVal();
  }

   if (chiSquared < bestChiSquared)
       bestChiSquared = chiSquared;
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

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
#include "OscLib/OscCalcPMNSOpt.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Systs/AnaSysts.h"
#include "CAFAna/Core/Loaders.h"

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TRandom3.h"

#include <random>

using namespace ana;

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/StateFilesDetectorSystematics.root";

void sensitivityFitSystematics(std::string &mode);

void PerformFit(MultiExperiment &experiment, std::vector<const ISyst*> &systematics, const double deltaCPSeed, const bool isLowerOctant, bool fitCPC, double &bestChiSquared);
std::vector<Spectrum> Get_Thrown_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);
std::vector<Spectrum> Get_DetectorSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, std::vector<std::string> &systematicsToShift, const float sigmaShift, const double deltaCP, const float pot);
Spectrum Get_CentralValue_NO(const PredictionInterp &predictionGenerator, const double deltaCP, const float pot);
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

std::default_random_engine generator;

int N_TEST_DELTA_CP_VALUES = 100; // 10
bool MAKE_THROWS = false;
int N_THROWS = MAKE_THROWS ? 200 : 1; // 500
bool FIT_SYSTEMATICS = false;

void sensitivityFitSystematics(int mode, float sigmaShift)
{

  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsSystematics_NO.root";
  vector<std::string> systematicsToShift;

  if (mode == 0)
  {
    std::cout << "All energy systematics" << std::endl;
    systematicsToShift = {
    "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt", 
    "ChargedHadUncorrFD", "UncorrFDHadSqrt", "UncorrFDHadInvSqrt", "ChargedHadResFD",
    "EScaleMuLArFD", "UncorrFDMuSqrt", "UncorrFDMuInvSqrt", "MuonResFD",
    "NUncorrFD", "UncorrFDNSqrt", "UncorrFDNInvSqrt", "NResFD",
    "EMUncorrFD", "UncorrFDEMSqrt", "UncorrFDEMInvSqrt", "EMResFD"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsAllEnergySystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else if (mode == 1)
  {
    std::cout << "Total energy scale and charged hadron systematics" << std::endl;

    systematicsToShift = {
    "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt", 
    "ChargedHadUncorrFD", "UncorrFDHadSqrt", "UncorrFDHadInvSqrt", "ChargedHadResFD"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsChargedHadronEnergySystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else if (mode == 2)
  {
    std::cout << "Total energy scale and muon systematics" << std::endl;

    systematicsToShift = {
    "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt", 
    "EScaleMuLArFD", "UncorrFDMuSqrt", "UncorrFDMuInvSqrt", "MuonResFD"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsMuonEnergySystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else if (mode == 3)
  {
    std::cout << "Total energy scale and neutron systematics" << std::endl;

    systematicsToShift = {
      "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt",
      "NUncorrFD", "UncorrFDNSqrt", "UncorrFDNInvSqrt", "NResFD"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsNeutronEnergySystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else if (mode == 4)
  {
    std::cout << "Total energy scale and EM systematics" << std::endl;

    systematicsToShift = {
      "EnergyScaleFD", "UncorrFDTotSqrt", "UncorrFDTotInvSqrt"
      "EMUncorrFD", "UncorrFDEMSqrt", "UncorrFDEMInvSqrt", "EMResFD"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsEMEnergySystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else if (mode == 5)
  {
    std::cout << "reconstruction model systematics" << std::endl;
    systematicsToShift = {"FDRecoNumuSyst", "FDRecoNueSyst"};

    outputFileName = "/storage/epp2/phrsnt/lblpwgtools/standardCAF/IZZLESensitivityPlotsRecoModelSystematics" + std::to_string(sigmaShift) + "SigmaShift_NO.root";
  }
  else
  {
    std::cout << "will apply no systematics shifts" << std::endl;
  }

  std::cout << "With " << sigmaShift << " sigma shift" << std::endl;


  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  PredictionInterp& interpGenNue_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_FHC_IZZLE").release();
  PredictionInterp& interpGenNue_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNue_RHC_IZZLE").release();
  PredictionInterp& interpGenNumu_FHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_FHC_IZZLE").release();
  PredictionInterp& interpGenNumu_RHC_IZZLE = *ana::LoadFrom<PredictionInterp>(inputFile, "interpGenNumu_RHC_IZZLE").release();

  std::vector<const PredictionInterp*> predictionGenerators = {&interpGenNue_FHC_IZZLE, &interpGenNue_RHC_IZZLE, &interpGenNumu_FHC_IZZLE, &interpGenNumu_RHC_IZZLE};

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

  for (int j = 0; j < nTestCPValues; ++j)
  {
      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);
      deltaCPValues.push_back(trueDeltaCP);

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j + 1) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i + 1) << "/" << std::to_string(nThrows) << std::endl;

          // Get the prediction
          std::vector<Spectrum> predictionVector;

	  if (MAKE_THROWS)
	  {
	    predictionVector = Get_Thrown_NO(predictionGenerators, trueDeltaCP, pot);
	  }
	  else
	  {
	    predictionVector = Get_DetectorSys_NO(predictionGenerators, systematicsToShift, sigmaShift, trueDeltaCP, pot);
	  }

          // Turn into an experiment (I tried to be fancy and it didn't work, leave me alone)
          MultiExperiment multiExperiment;
          const SingleSampleExperiment experimentNueFHC(predictionGenerators[0], predictionVector[0]);
          const SingleSampleExperiment experimentNueRHC(predictionGenerators[1], predictionVector[1]);
          const SingleSampleExperiment experimentNumuFHC(predictionGenerators[2], predictionVector[2]);
          const SingleSampleExperiment experimentNumuRHC(predictionGenerators[3], predictionVector[3]);
          multiExperiment.Add(&experimentNueFHC);
          multiExperiment.Add(&experimentNueRHC);
          multiExperiment.Add(&experimentNumuFHC);
          multiExperiment.Add(&experimentNumuRHC);

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
          std::cout << "significance: " << TMath::Sqrt(std::fabs(chiSquaredDifference)) << std::endl;
          std::cout << "////////////////////////////" << std::endl;

	  chiSquaredCPCValues.push_back(bestChiSquaredCPC);
	  chiSquaredCPVValues.push_back(bestChiSquaredCPV);
          chiSquaredValues.push_back(chiSquaredDifference);
      }

      tree->Branch("chiSquaredCPCValues" , &chiSquaredCPCValues);
      tree->Branch("chiSquaredCPVValues", &chiSquaredCPVValues);
      tree->Branch("chiSquaredValues", &chiSquaredValues);
      tree->Branch("deltaCPValues", &deltaCPValues);  
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

std::vector<Spectrum> Get_DetectorSys_NO(std::vector<const PredictionInterp*> &predictionGenerators, std::vector<std::string> &systematicsToShift, const float sigmaShift, const double deltaCP, const float pot)
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
	      std::cout << "yoy yo" << std::endl;
              std::cout << "ShortName: " << systematic->ShortName() << std::endl;
	      systematicShiftMap[systematic] = sigmaShift;
	  }
      }
    }

    const SystShifts systematicShifts(systematicShiftMap);

    // Make oscillation throw
    osc::IOscCalcAdjustable* calc = DefaultOscCalc();
    SetOscCalc_CentralValue_NO(calc, deltaCP);

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

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

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

#include <random>

// also add in some RHC files!

using namespace ana;

const std::string INPUT_FILE_NAME = "/dune/data/users/imawby/standardCAF/StateFiles.root";

void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment, const double deltaCPSeed, const bool isLowerOctant, double &bestChiSquared, bool fitCPC);
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

std::default_random_engine generator;

bool MAKE_THROWS = false;
bool N_THROWS = MAKE_THROWS ? 200 : 1;

void sensitivityFit()
{
  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  PredictionNoExtrap& predNue_FHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNue_FHC_IZZLE").release();
  PredictionNoExtrap& predNumu_FHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNumu_FHC_IZZLE").release();
  PredictionNoExtrap& predNue_RHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNue_RHC_IZZLE").release();
  PredictionNoExtrap& predNumu_RHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNumu_RHC_IZZLE").release();

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile("/dune/data/users/imawby/standardCAF/IZZLESensitivityPlots_NO.root", "CREATE");

  std::cout << "created output file" << std::endl;

  //int nTestCPValues(10);
  //unsigned int nThrows(500);

  int nTestCPValues(2);
  unsigned int nThrows(2);

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

          // Make oscillation throw 
          osc::IOscCalcAdjustable* calc = DefaultOscCalc();
          MAKE_THROWS ? SetOscCalc_Throw_NO(calc, trueDeltaCP) : SetOscCalc_CentralValue_NO(calc, trueDeltaCP);

          // Make the mock data for trueDeltaCP value
          // No statistical fluctuations
          const Spectrum nueFluctuations_FHC = predNue_FHC_IZZLE.Predict(calc).AsimovData(pot);
          const Spectrum numuFluctuations_FHC = predNumu_FHC_IZZLE.Predict(calc).AsimovData(pot);
          const Spectrum nueFluctuations_RHC = predNue_RHC_IZZLE.Predict(calc).AsimovData(pot);
          const Spectrum numuFluctuations_RHC = predNumu_RHC_IZZLE.Predict(calc).AsimovData(pot);

          // Apply statistical fluctuations
          /*
          const Spectrum nueFluctuations_FHC = predNue_FHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum numuFluctuations_FHC = predNumu_FHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum nueFluctuations_RHC = predNue_RHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum numuFluctuations_RHC = predNumu_RHC_IZZLE.Predict(calc).MockData(pot);
          */

          // Create an experiment object to compare predictions to my data
          SingleSampleExperiment nueExperiment_FHC(&predNue_FHC_IZZLE, nueFluctuations_FHC);
          SingleSampleExperiment numuExperiment_FHC(&predNumu_FHC_IZZLE, numuFluctuations_FHC);
          SingleSampleExperiment nueExperiment_RHC(&predNue_RHC_IZZLE, nueFluctuations_RHC);
          SingleSampleExperiment numuExperiment_RHC(&predNumu_RHC_IZZLE, numuFluctuations_RHC);
          //MultiExperiment multiExperiment({&nueExperiment_FHC, &numuExperiment_FHC, &nueExperiment_RHC, &numuExperiment_RHC});
          MultiExperiment multiExperiment({&nueExperiment_FHC, &numuExperiment_FHC});

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          osc::IOscCalcAdjustable* fitCalc = DefaultOscCalc();
          double bestChiSquaredCPC(std::numeric_limits<float>::max());
          double bestChiSquaredCPV(std::numeric_limits<float>::max());

          // CPC Fits
          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), true, bestChiSquaredCPC, true);
          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), false, bestChiSquaredCPC, true);

          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), true, bestChiSquaredCPC, true);
          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), false, bestChiSquaredCPC, true);

          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), true, bestChiSquaredCPC, true);
          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), false, bestChiSquaredCPC, true);

          // CPV Fits
          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), true, bestChiSquaredCPV, false);
          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), false, bestChiSquaredCPV, false);

          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), true, bestChiSquaredCPV, false);
          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), false, bestChiSquaredCPV, false);

          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), true, bestChiSquaredCPV, false);
          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), false, bestChiSquaredCPV, false);

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

void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment, const double deltaCPSeed, const bool isLowerOctant, double &bestChiSquared, bool fitCPC)
{
  SetOscCalc_CentralValue_NO(calc, deltaCPSeed);

  double chiSquared(9999);

  std::vector<const IFitVar *> fitVariables = {&kFitDmSq32NHScaled};

  isLowerOctant ? fitVariables.push_back(&kFitSinSqTheta23LowerOctant) : fitVariables.push_back(&kFitSinSqTheta23UpperOctant);

  if (!fitCPC)
      fitVariables.push_back(&kFitDeltaInPiUnits);

  if (isLowerOctant)
  {
      calc->SetTh23(0.735);
      MinuitFitter fit(&experiment, fitVariables);
      chiSquared = fit.Fit(calc)->EvalMetricVal();
  }
  else
  {
      calc->SetTh23(0.835);
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

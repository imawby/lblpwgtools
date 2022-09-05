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

void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment,
  const double deltaCPSeed, const bool isLowerOctant, double &bestDeltaCP, double &bestChiSquared);
void SetOscCalc_CentralValue_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_NO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_CentralValue_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void SetOscCalc_Throw_IO(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

std::default_random_engine generator;

void deltaCPResolution_statisticalOnly()
{
  std::cout << "Reading caf files..." << std::endl;

  TFile * inputFile = TFile::Open(INPUT_FILE_NAME.c_str());

  PredictionNoExtrap& predNue_FHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNue_FHC_IZZLE").release();
  PredictionNoExtrap& predNumu_FHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNumu_FHC_IZZLE").release();
  PredictionNoExtrap& predNue_RHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNue_RHC_IZZLE").release();
  PredictionNoExtrap& predNumu_RHC_IZZLE = *ana::LoadFrom<PredictionNoExtrap>(inputFile, "predNumu_RHC_IZZLE").release();

  inputFile->Close();

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  TFile * outputFile = new TFile("/dune/data/users/imawby/standardCAF/IZZLEResolutionPlots_NO.root", "CREATE");

  std::cout << "created output file" << std::endl;

  int nTestCPValues(10);
  unsigned int nThrows(500);
  float stepSizeCP((2.0 * TMath::Pi()) / static_cast<float>(nTestCPValues));

  double meanResolutionValues[nTestCPValues], oneSigmaResolutionValues[nTestCPValues], twoSigmaResolutionValues[nTestCPValues], threeSigmaResolutionValues[nTestCPValues];
  double testCPValues[nTestCPValues];

  for (int j = 0; j < nTestCPValues; ++j)
  {
      const double trueDeltaCP(static_cast<float>(j) * stepSizeCP);

      TH1D * bestFitDeltaCPValues = new TH1D(("bestFitDeltaCPValues_" + std::to_string(trueDeltaCP)).c_str(), ("bestFitDeltaCPValues_" + std::to_string(trueDeltaCP)).c_str(), 
          41, -0.025 * TMath::Pi(), 2.025 * TMath::Pi());

      for (unsigned int i = 0; i < nThrows; ++i)
      {
          std::cout << "iteration: " << std::to_string(j) << "/" << std::to_string(nTestCPValues) << std::endl;
          std::cout << "nThrow: " << std::to_string(i) << "/" << std::to_string(nThrows) << std::endl;

          // Make oscillation throw 
          osc::IOscCalcAdjustable* calc = DefaultOscCalc();
          SetOscCalc_Throw_NO(calc, trueDeltaCP);

          // Make the mock data for trueDeltaCP value
          const Spectrum nueFluctuations_FHC = predNue_FHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum numuFluctuations_FHC = predNumu_FHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum nueFluctuations_RHC = predNue_RHC_IZZLE.Predict(calc).MockData(pot);
          const Spectrum numuFluctuations_RHC = predNumu_RHC_IZZLE.Predict(calc).MockData(pot);

          // Create an experiment object to compare predictions to my data
          SingleSampleExperiment nueExperiment_FHC(&predNue_FHC_IZZLE, nueFluctuations_FHC);
          SingleSampleExperiment numuExperiment_FHC(&predNumu_FHC_IZZLE, numuFluctuations_FHC);
          SingleSampleExperiment nueExperiment_RHC(&predNue_RHC_IZZLE, nueFluctuations_RHC);
          SingleSampleExperiment numuExperiment_RHC(&predNumu_RHC_IZZLE, numuFluctuations_RHC);
          MultiExperiment multiExperiment({&nueExperiment_FHC, &numuExperiment_FHC, &nueExperiment_RHC, &numuExperiment_RHC});

          // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
          osc::IOscCalcAdjustable* fitCalc = DefaultOscCalc();
          double bestChiSquared(std::numeric_limits<float>::max()), bestDeltaCP(999);

          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), true, bestDeltaCP, bestChiSquared);
          PerformFit(fitCalc, multiExperiment, 0.0 * TMath::Pi(), false, bestDeltaCP, bestChiSquared);

          PerformFit(fitCalc, multiExperiment, 0.5 * TMath::Pi(), true, bestDeltaCP, bestChiSquared);
          PerformFit(fitCalc, multiExperiment, 0.5 * TMath::Pi(), false, bestDeltaCP, bestChiSquared);

          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), true, bestDeltaCP, bestChiSquared);
          PerformFit(fitCalc, multiExperiment, 1.0 * TMath::Pi(), false, bestDeltaCP, bestChiSquared);

          PerformFit(fitCalc, multiExperiment, 1.5 * TMath::Pi(), true, bestDeltaCP, bestChiSquared);
          PerformFit(fitCalc, multiExperiment, 1.5 * TMath::Pi(), false, bestDeltaCP, bestChiSquared);

          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), true, bestDeltaCP, bestChiSquared);
          PerformFit(fitCalc, multiExperiment, 2.0 * TMath::Pi(), false, bestDeltaCP, bestChiSquared);

          bestFitDeltaCPValues->Fill(bestDeltaCP);

          /*
          std::cout << "///////////////////////////////" << std::endl;
          std::cout << "i: " << i << std::endl;
          std::cout << "True DeltaCP: " << trueDeltaCP << std::endl;
          std::cout << "Best DeltaCP: " << bestDeltaCP << std::endl;
          std::cout << "///////////////////////////////" << std::endl;
          */
      }

      // Now need to shift histogram so i can fit a gaussian to a periodic function

      double bestEntries(-9999);
      double bestCentre(9999);

      for (int i = 1; i <= bestFitDeltaCPValues->GetNbinsX(); ++i)
      {
          const double binContent(bestFitDeltaCPValues->GetBinContent(i));
          if (binContent > bestEntries)
          {
              bestEntries = binContent;
              bestCentre = bestFitDeltaCPValues->GetBinCenter(i);
          }
      }

      const double min(bestCentre - (1.025 * TMath::Pi()));
      const double max(bestCentre + (1.025 * TMath::Pi()));

      TH1D * shiftedDeltaCPValues = new TH1D(("shiftedDeltaCPValues_" + std::to_string(trueDeltaCP)).c_str(), ("shiftedDeltaCPValues_" + std::to_string(trueDeltaCP)).c_str(), 
          41, min, max);

      for (int i = 1; i <= bestFitDeltaCPValues->GetNbinsX(); ++i)
      {
          const float binCenter(bestFitDeltaCPValues->GetBinCenter(i));
          const float binContent(bestFitDeltaCPValues->GetBinContent(i));

          if (binContent < std::numeric_limits<float>::epsilon())
              continue;

          const double separation(binCenter - bestCentre);
          int shiftedBin(GetBin(shiftedDeltaCPValues, binCenter));

          if (std::fabs(separation) > TMath::Pi())
          {
              if (separation > 0.0)
              {
                  double jam = binCenter - (2.0 * TMath::Pi());
                  shiftedBin = GetBin(shiftedDeltaCPValues, jam);
              }
              else
              {
                  double jam = binCenter + (2.0 * TMath::Pi());
                  shiftedBin = GetBin(shiftedDeltaCPValues, jam);
              }
          }

          shiftedDeltaCPValues->SetBinContent(shiftedBin, binContent);
      }

      shiftedDeltaCPValues->Scale(1.0 / shiftedDeltaCPValues->Integral());

      // now fit a gaussian!!! 
      TF1 * gaussianFit = new TF1(("gaussianFit_" + std::to_string(j)).c_str() , "(1.0/([0]*sqrt(2.0*3.14))) * exp(-0.5 * ((x - [1])/ [0])^2)", min, max);
      gaussianFit->SetParameter(1, bestCentre);
      shiftedDeltaCPValues->Fit(("gaussianFit_" + std::to_string(j)).c_str());

      float mean = gaussianFit->GetParameter(1);
      float standardDeviation = gaussianFit->GetParameter(0);

      bestFitDeltaCPValues->Write(("bestFitDeltaCPValues_" + std::to_string(j)).c_str());
      shiftedDeltaCPValues->Write(("shiftedDeltaCPValues_" + std::to_string(j)).c_str());
      gaussianFit->Write(("gaussianFit" + std::to_string(j)).c_str());

      testCPValues[j] = (trueDeltaCP * 180.0) / TMath::Pi();

      float resolution = std::min(std::fabs(mean - trueDeltaCP), std::fabs(mean - (trueDeltaCP + (2.0 * TMath::Pi()))));

      // convert to degrees because i am dumb

      standardDeviation *= (180.0 / TMath::Pi());
      resolution *= (180.0 / TMath::Pi());

      meanResolutionValues[j] = resolution;
      oneSigmaResolutionValues[j] = resolution + (standardDeviation);
      twoSigmaResolutionValues[j] = resolution + (2.0 * standardDeviation);
      threeSigmaResolutionValues[j] = resolution + (3.0 * standardDeviation);
  }

  TGraph * meanResolutionGraph = new TGraph(nTestCPValues, testCPValues, meanResolutionValues);
  TGraph * oneSigmaResolutionGraph = new TGraph(nTestCPValues, testCPValues, oneSigmaResolutionValues);
  TGraph * twoSigmaResolutionGraph = new TGraph(nTestCPValues, testCPValues, twoSigmaResolutionValues);
  TGraph * threeSigmaResolutionGraph = new TGraph(nTestCPValues, testCPValues, threeSigmaResolutionValues);

  meanResolutionGraph->Write("meanResolutionGraph");
  oneSigmaResolutionGraph->Write("oneSigmaResolutionGraph");
  twoSigmaResolutionGraph->Write("twoSigmaResolutionGraph");
  threeSigmaResolutionGraph->Write("threeSigmaResolutionGraph");
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


void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment,
  const double deltaCPSeed, const bool isLowerOctant, double &bestDeltaCP, double &bestChiSquared)
{
  SetOscCalc_CentralValue_NO(calc, deltaCPSeed);

  double deltaCP(9999), chiSquared(9999);

  if (isLowerOctant)
  {
      calc->SetTh23(0.735);
      MinuitFitter fit(&experiment, {&kFitDmSq32NHScaled, &kFitSinSqTheta23LowerOctant, &kFitDeltaInPiUnits});
      chiSquared = fit.Fit(calc)->EvalMetricVal();
      deltaCP = calc->GetdCP();
  }
  else
  {
      calc->SetTh23(0.835);
      MinuitFitter fit(&experiment, {&kFitDmSq32NHScaled, &kFitSinSqTheta23UpperOctant, &kFitDeltaInPiUnits});
      chiSquared = fit.Fit(calc)->EvalMetricVal();
      deltaCP = calc->GetdCP();
  }

   if (chiSquared < bestChiSquared)
   {
       bestChiSquared = chiSquared;
       bestDeltaCP = deltaCP;
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

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

// also add in some RHC files!

using namespace ana;

void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment,
  const double deltaCPSeed, const bool isLowerOctant, double &bestDeltaCP, double &bestChiSquared);
void SetOscCalc_NuFit(osc::IOscCalcAdjustable *& calc, const double deltaCP);
void deltaCPResolution_statisticalOnly();
int GetBin(TH1D * hist, double value);

void deltaCPResolution_statisticalOnly()
{
  std::cout << "Reading caf files..." << std::endl;

  const std::string fnameNonSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nu.root";
  const std::string fnameNueSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nue.root";
  const std::string fnameTauSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nutau.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);

  std::cout << "Read files!" << std::endl;

  const Var kRecoNueEnergy = SIMPLEVAR(Ev_reco_nue);
  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);

  const Var kTrackPandizzle = SIMPLEVAR(selTrackPandizzleScore);
  const Var kShowerPandrizzle = SIMPLEVAR(selShowerPandrizzleScore);

  const Cut kPassesCVN_nue = SIMPLEVAR(cvnnue) > 0.85;
  const Cut kPassesCVN_numu = SIMPLEVAR(cvnnumu) > 0.5;

  // TODO - change the binning
  const Binning binsEnergy = Binning::Simple(80, 0, 10);
  const HistAxis axNueEnergy("Reco #nu_{e} energy (GeV)", binsEnergy, kRecoNueEnergy);
  const HistAxis axNumuEnergy("Reco #nu_{#mu} energy (GeV)", binsEnergy, kRecoNumuEnergy);

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap predNue(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kIsNueSelected);
  PredictionNoExtrap predNumu(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kIsNumuSelected);
  //PredictionNoExtrap predNue(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNueEnergy, kPassesCVN_nue);
  //PredictionNoExtrap predNumu(loaderNonSwap, loaderNueSwap, loaderTauSwap, axNumuEnergy, kPassesCVN_numu);

  std::cout << "Filling spectra..." << std::endl;

  // These calls will fill all of the constituent parts of the prediction
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();

  std::cout << "Spectra filled!" << std::endl;

  const double trueDeltaCP(0.0 * TMath::Pi());
  TH1D * bestFitDeltaCPValues = new TH1D("bestFitDeltaCPValues", "bestFitDeltaCPValues", 100, 0.0, 2.0 * M_PI);
  bestFitDeltaCPValues->SetDirectory(0);

  // LOOP OVER DELTA CP VALUES
      osc::IOscCalcAdjustable* calc = DefaultOscCalc();
      SetOscCalc_NuFit(calc, trueDeltaCP);

      // Make the mock data for trueDeltaCP value
      const Spectrum nueFluctuations_FHC = predNue.Predict(calc).MockData(pot);
      const Spectrum numuFluctuations_FHC = predNumu.Predict(calc).MockData(pot);

      // Create an experiment object to compare predictions to my data
      SingleSampleExperiment nueExperiment_FHC(&predNue, nueFluctuations_FHC);
      SingleSampleExperiment numuExperiment_FHC(&predNumu, numuFluctuations_FHC);
      MultiExperiment multiExperiment({&nueExperiment_FHC, &numuExperiment_FHC});

      // Make several fits to avoid falling into the wrong minima (find the best deltaCP)
      double bestChiSquared(std::numeric_limits<float>::max()), bestDeltaCP(999);

      PerformFit(calc, multiExperiment, -1.0 * M_PI, true, bestDeltaCP, bestChiSquared);
      PerformFit(calc, multiExperiment, -1.0 * M_PI, false, bestDeltaCP, bestChiSquared);

      PerformFit(calc, multiExperiment, -0.5 * M_PI, true, bestDeltaCP, bestChiSquared);
      PerformFit(calc, multiExperiment, -0.5 * M_PI, false, bestDeltaCP, bestChiSquared);

      PerformFit(calc, multiExperiment, 0.0 * M_PI, true, bestDeltaCP, bestChiSquared);
      PerformFit(calc, multiExperiment, 0.0 * M_PI, false, bestDeltaCP, bestChiSquared);

      PerformFit(calc, multiExperiment, 0.5 * M_PI, true, bestDeltaCP, bestChiSquared);
      PerformFit(calc, multiExperiment, 0.5 * M_PI, false, bestDeltaCP, bestChiSquared);

      bestDeltaCP -= TMath::Pi();

      SetOscCalc_NuFit(calc, trueDeltaCP);
      MinuitFitter fit(&multiExperiment, {&kFitDmSq32NHScaled, &kFitSinSqTheta23});

      Resolution * resolution = new Resolution(&fit, calc, 0, bestChiSquared);

      double functionMin = -1.5 * TMath::Pi();
      double functionMax = 1.5 * TMath::Pi();

      TF1 * fitFunction = new TF1("f", resolution, &Resolution::FitResult, functionMin, functionMax, 0, "Resolution", "FitResult"); 

      // Search for two minima in the [true deltaCP, true deltaCP + pi] and [true deltaCP - pi, true deltaCP] ranges
      ROOT::Math::WrappedTF1 wrappedFitFunction(*fitFunction);
      const double precision1 = 1e-5, precision2 = 1e-7;

      ROOT::Math::BrentRootFinder upperRootFinder;
      upperRootFinder.SetFunction(wrappedFitFunction, bestDeltaCP, bestDeltaCP + (TMath::Pi() / 2.0)); // Set function and interval in which to search for root
      upperRootFinder.SetNpx(15);                                                                      // Set the number of point used to bracket root using a grid
      upperRootFinder.Solve(100, precision1, precision2);                                              // The grid is searched to find the minimum position
      double upperRoot = upperRootFinder.Root();
      double upperSeparation = TMath::Abs(upperRoot - trueDeltaCP);

      ROOT::Math::BrentRootFinder lowerRootFinder;
      lowerRootFinder.SetFunction(wrappedFitFunction, bestDeltaCP - (TMath::Pi() / 2.0), bestDeltaCP);
      lowerRootFinder.SetNpx(15);
      lowerRootFinder.Solve(100, precision1, precision2);
      double lowerRoot = lowerRootFinder.Root();
      double lowerSeparation = TMath::Abs(lowerRoot - trueDeltaCP);

      //double averageResolution = (upperSeparation + lowerSeparation) / 2.0; // work out an average?
      
      std::cout << "///////////////////////////////" << std::endl;
      std::cout << "i: " << i << std::endl;
      std::cout << "True DeltaCP: " << trueDeltaCP << std::endl;
      std::cout << "Lower Fit DeltaCP: " << lowerRoot << std::endl;
      std::cout << "Lower Resolution: " << lowerSeparation << std::endl;
      std::cout << "Upper Fit DeltaCP: " << upperRoot << std::endl;
      std::cout << "Upper Resolution: " << upperSeparation << std::endl;
      //std::cout << "Average Resolution: " << averageResolution << std::endl;
      std::cout << "///////////////////////////////" << std::endl;

      //maybe not?? this just gives one resolution value?


  TCanvas * bee = new TCanvas("bee", "bee");
  bee->Divide(2,1);
  bee->cd(1);
  bestFitDeltaCPValues->SetTitle("[0, 2pi];Resolution;nEntries");
  bestFitDeltaCPValues->SetLineColor(kPink);
  bestFitDeltaCPValues->SetLineWidth(2);
  bestFitDeltaCPValues->Draw("hist");
  bee->cd(2);
  shiftedDeltaCPValues->SetTitle("shifted;Resolution;nEntries");
  shiftedDeltaCPValues->SetLineColor(kPink);
  shiftedDeltaCPValues->SetLineWidth(2);
  shiftedDeltaCPValues->GetXaxis()->SetRangeUser(min, max);
  shiftedDeltaCPValues->Draw("hist");
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

void SetOscCalc_NuFit(osc::IOscCalcAdjustable *& calc, const double deltaCP)
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

void PerformFit(osc::IOscCalcAdjustable *& calc, MultiExperiment &experiment,
  const double deltaCPSeed, const bool isLowerOctant, double &bestDeltaCP, double &bestChiSquared)
{
  SetOscCalc_NuFit(calc, deltaCPSeed);

  double deltaCP(9999), chiSquared(9999);

  if (isLowerOctant)
  {
      MinuitFitter fit(&experiment, {&kFitDmSq32NHScaled, &kFitSinSqTheta23LowerOctant, &kFitDeltaInPiUnits});
      chiSquared = fit.Fit(calc)->EvalMetricVal();
      deltaCP = calc->GetdCP();
  }
  else
  {
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

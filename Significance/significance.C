// Make significance plot
// cafe significance.C

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

using namespace ana;

void significance()
{

  std::cout << "Hello Izzie!" << std::endl;

  const std::string fnameNonSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nu.root";
  const std::string fnameNueSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nue.root";
  const std::string fnameTauSwap = "/dune/app/users/imawby/selection/cafs/caf_tree_nutau.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);

  const Var kRecoNueEnergy = SIMPLEVAR(Ev_reco_nue);
  const Var kRecoNumuEnergy = SIMPLEVAR(Ev_reco_numu);

  const Var kTrackPandizzle = SIMPLEVAR(selTrackPandizzleScore);
  const Var kShowerPandrizzle = SIMPLEVAR(selShowerPandrizzleScore);

  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axNueEnergy("Reco #nu_{e} energy (GeV)", binsEnergy, kRecoNueEnergy);
  const HistAxis axNumuEnergy("Reco #nu_{#mu} energy (GeV)", binsEnergy, kRecoNumuEnergy);

  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap predNue(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                             axNueEnergy, kIsNueSelected);

  PredictionNoExtrap predNumu(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                              axNumuEnergy, kIsNumuSelected);

  // These calls will fill all of the constituent parts of the prediction
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();

  osc::IOscCalcAdjustable* calc = DefaultOscCalc();

  calc->SetdCP(0);
  const Spectrum zeroSpectrumNue = predNue.Predict(calc).AsimovData(pot);
  const Spectrum zeroSpectrumNumu = predNumu.Predict(calc).AsimovData(pot);

  calc->SetdCP(TMath::Pi());
  const Spectrum piSpectrumNue = predNue.Predict(calc).AsimovData(pot);
  const Spectrum piSpectrumNumu = predNumu.Predict(calc).AsimovData(pot);

  const unsigned int nSteps = 50;
  const float stepSize = (2.0 * TMath::Pi()) / static_cast<float>(nSteps);

  Double_t deltaCP[nSteps];
  Double_t sigNue[nSteps], sigNumu[nSteps], sigTotal[nSteps];

  for (unsigned int i = 0; i < nSteps; ++i)
  {
      const float deltaCPValue = ((-1) * TMath::Pi()) + (i * stepSize);
      std::cout << "deltaCP: " << deltaCPValue << std::endl;
      calc->SetdCP(deltaCPValue);

      SingleSampleExperiment zeroExpNue(&predNue, zeroSpectrumNue);
      const float zeroChiSqNue = zeroExpNue.ChiSq(calc); 

      SingleSampleExperiment zeroExpNumu(&predNumu, zeroSpectrumNumu);
      const float zeroChiSqNumu = zeroExpNumu.ChiSq(calc); 

      if (std::fabs(deltaCPValue - calc->GetdCP()) > std::numeric_limits<float>::epsilon())
          std::cout << "DELTA CP HAS CHANGED" << std::endl;

      SingleSampleExperiment piExpNue(&predNue, piSpectrumNue);
      const float piChiSqNue = piExpNue.ChiSq(calc);

      SingleSampleExperiment piExpNumu(&predNumu, piSpectrumNumu);
      const float piChiSqNumu = piExpNumu.ChiSq(calc);

      if (std::fabs(deltaCPValue - calc->GetdCP()) > std::numeric_limits<float>::epsilon())
          std::cout << "DELTA CP HAS CHANGED" << std::endl;

      double sigNueValue = std::sqrt(std::min(zeroChiSqNue, piChiSqNue));
      double sigNumuValue = std::sqrt(std::min(zeroChiSqNumu, piChiSqNumu));
      double sigTotalValue = std::sqrt(std::min(zeroChiSqNue, piChiSqNue) + std::min(zeroChiSqNumu, piChiSqNumu));



      deltaCP[i] = deltaCPValue;

      std::cout << "deltaCP: " << deltaCP[i] << ", sigNue: " << sigNueValue << std::endl;
      sigNue[i] = sigNueValue;
      sigNumu[i] = sigNumuValue;
      sigTotal[i] = sigTotalValue;
  }

  std::cout << "////////////////////////" << std::endl;
  std::cout << "Oscillation Calc Params: " << std::endl;
  std::cout << "////////////////////////" << std::endl;
  std::cout << "Baseline: " << calc->GetL() << std::endl;
  std::cout << "Rho: " << calc->GetRho() << std::endl;
  std::cout << "Dmsq21: " << calc->GetDmsq21() << std::endl;
  std::cout << "Dmsq32: " << calc->GetDmsq32() << std::endl;
  std::cout << "Theta12: " << calc->GetTh12() << std::endl;
  std::cout << "Theta13: " << calc->GetTh13() << std::endl;
  std::cout << "Theta23: " << calc->GetTh23() << std::endl;
  std::cout << "DeltaCP: " << calc->GetdCP() << std::endl; 

  TCanvas * jam = new TCanvas("jam", "jam");
  zeroSpectrumNumu.ToTH1(pot, kBlue)->Draw("hist same");

  TCanvas * cSigNue = new TCanvas("cSigNue", "cSigNue", 200, 10, 500, 300);
  TGraph* grNue = new TGraph(nSteps,deltaCP,sigNue);
  grNue->SetTitle("Nue;deltaCP [radians]; #sqrt{#chi^{2}}");;
  grNue->Draw("AC*");

  TCanvas * cSigNumu = new TCanvas("cSigNumu", "cSigNumu", 200, 10, 500, 300);
  TGraph* grNumu = new TGraph(nSteps,deltaCP,sigNumu);
  grNumu->SetTitle("Numu;deltaCP [radians]; #sqrt{#chi^{2}}");;
  grNumu->Draw("AC*");

  TCanvas * cSigTotal = new TCanvas("cSigTotal", "cSigTotal", 200, 10, 500, 300);
  TGraph* grTotal = new TGraph(nSteps,deltaCP,sigTotal);
  grTotal->SetTitle("Total;deltaCP [radians]; #sqrt{#chi^{2}}");;
  grTotal->Draw("AC*");
}

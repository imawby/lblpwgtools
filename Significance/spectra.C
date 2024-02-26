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

const std::string INPUT_FILE_NAME = "/storage/epp2/phrsnt/lblpwgtools/realRecoStandardCAF/fullEstimate/StateFilesAllSystematicsSplitBySign.root";

void spectra();

std::vector<Spectrum> Get_Spectra_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot);

void spectra()
{
  std::string outputFileName = "/storage/epp2/phrsnt/lblpwgtools/realRecoStandardCAF/spectra.root";

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

  const double pot = 3.5 * 1.1e21 * 40/1.13;

  TFile * outputFile = new TFile(outputFileName.c_str(), "CREATE");
  std::cout << "created output file" << std::endl;

  for (int j = 0; j < 11; ++j)
  {
    const double trueDeltaCP(static_cast<float>(j) * 0.1 * TMath::Pi());

    // Get the prediction
    std::vector<Spectrum> predictionVector;
    predictionVector = Get_Spectra_NO(predictionGenerators, trueDeltaCP, pot);

    std::string nueName = "nue_" + std::to_string(j);
    std::string anueName = "anue_" + std::to_string(j);
    std::string numuName = "numu_" + std::to_string(j);
    std::string anumuName = "anumu_" + std::to_string(j);

    predictionVector[0].ToTH1(pot)->Write(nueName.c_str());
    predictionVector[1].ToTH1(pot)->Write(anueName.c_str());
    predictionVector[2].ToTH1(pot)->Write(numuName.c_str());
    predictionVector[3].ToTH1(pot)->Write(anumuName.c_str());
  }

  outputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<Spectrum> Get_Spectra_NO(std::vector<const PredictionInterp*> &predictionGenerators, const double deltaCP, const float pot)
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


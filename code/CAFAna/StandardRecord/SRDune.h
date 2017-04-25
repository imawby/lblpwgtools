#pragma once

namespace caf
{
  class SRDune
  {
  public:
    // Reco info
    double Ev_reco;
    double mvaresult;

    // Truth info
    double Ev;
    //  float enu_truth; // so what's this one?
    int ccnc;
    int beamPdg;
    int neu;

    int run;

    int lep;
    int scat;
    int nipiz;
    int nipip;
    int nipim;
    float Q2;
  };
}

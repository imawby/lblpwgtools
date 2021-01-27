#pragma once

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/Var.h"

#include "fhiclcpp/ParameterSet.h"

#include <array>
#include <iostream>
#include <string>
#include <utility>

namespace PRISM {
enum class NuChan;
enum class BeamMode;
struct BeamChan;
struct MatchChan;

bool operator&(NuChan const &l, NuChan const &r);
std::ostream &operator<<(std::ostream &os, NuChan const &nc);
std::ostream &operator<<(std::ostream &os, BeamMode const &bm);
std::ostream &operator<<(std::ostream &os, BeamChan const &bm);
bool operator==(BeamChan const &l, BeamChan const &r);
bool operator<(BeamChan const &l, BeamChan const &r);
bool operator<(MatchChan const &l, MatchChan const &r);
std::ostream &operator<<(std::ostream &os, MatchChan const &mc);

enum class NuChan {
  kNumuApp = (1 << 0),
  kNumuBarApp = (1 << 1),
  kNueApp = (1 << 2),
  kNueBarApp = (1 << 3),

  kNumuIntrinsic = (1 << 4),
  kNumuBarIntrinsic = (1 << 5),
  kNueIntrinsic = (1 << 6),
  kNueBarIntrinsic = (1 << 7),

  kNumu = kNumuApp | kNumuIntrinsic,
  kNumuBar = kNumuBarApp | kNumuBarIntrinsic,
  kNue = kNueApp | kNueIntrinsic,
  kNueBar = kNueBarApp | kNueBarIntrinsic,

  kNumuNumuBar = kNumu | kNumuBar,
  kNueNueBar = kNue | kNueBar
};

enum class BeamMode {
  kNuMode_293kA = (1 << 0),
  kNuMode_280kA = (1 << 1),
  kNuMode = kNuMode_293kA | kNuMode_280kA,
  kNuBarMode_293kA = (1 << 2),
  kNuBarMode_280kA = (1 << 3),
  kNuBarMode = kNuBarMode_293kA | kNuBarMode_280kA
};

struct BeamChan {
  BeamMode mode;
  NuChan chan;
};

static BeamChan const kNumu_Numode = {BeamMode::kNuMode,
                                      NuChan::kNumuIntrinsic};
static BeamChan const kNumuBar_NuBarmode = {BeamMode::kNuBarMode,
                                            NuChan::kNumuBarIntrinsic};
static BeamChan const kNue_Numode = {BeamMode::kNuMode, NuChan::kNueApp};
static BeamChan const kNueBar_NuBarmode = {BeamMode::kNuBarMode,
                                           NuChan::kNueBarApp};

inline BeamChan GetBeamChan(std::string const &descript, bool IsND) {

  if (IsND) {
    if (descript == "numu_numode") {
      return kNumu_Numode;
    } else if (descript == "numubar_nubarmode") {
      return kNumuBar_NuBarmode;
    } else {
      std::cout << "[ERROR]: Invalid ND beam-chan, currently support "
                   "\"numu_numode\" or \"numubar_nubarmode\"."
                << std::endl;
      abort();
    }
  } else {
    if (descript == "numu_numode") {
      return kNumu_Numode;
    } else if (descript == "numubar_nubarmode") {
      return kNumuBar_NuBarmode;
    } else if (descript == "nue_numode") {
      return kNue_Numode;
    } else if (descript == "nuebar_nubarmode") {
      return kNueBar_NuBarmode;
    } else {
      std::cout << "[ERROR]: Invalid ND beam-chan, currently support "
                   "\"numu_numode\", \"numubar_nubarmode\", \"nue_numode\", "
                   "\"nuebar_nubarmode\"."
                << std::endl;
      abort();
    }
  }
}

struct MatchChan {
  BeamChan from;
  BeamChan to;
};

static MatchChan const kNumuDisappearance_Numode = {kNumu_Numode, kNumu_Numode};
static MatchChan const kNumuBarDisappearance_NuBarmode = {kNumuBar_NuBarmode,
                                                          kNumuBar_NuBarmode};
static MatchChan const kNueAppearance_Numode = {kNumu_Numode, kNue_Numode};
static MatchChan const kNueBarAppearance_NuBarmode = {kNumuBar_NuBarmode,
                                                      kNueBar_NuBarmode};

inline std::string GetMatchChanShortName(MatchChan ch) {
  if (ch.to.mode == BeamMode::kNuMode) {
    return (ch.to.chan & NuChan::kNumu) ? "NumuDisp" : "NueApp";
  } else {
    return (ch.to.chan & NuChan::kNumuBar) ? "NumuBarDisp" : "NueBarApp";
  }
  return "UnknownChannel";
}

inline MatchChan GetMatchChan(fhicl::ParameterSet const &ps) {
  return {GetBeamChan(ps.get<std::string>("ND"), true),
          GetBeamChan(ps.get<std::string>("FD"), false)};
}

// Enum-like list of Ids for use in lists of PRISM objects
size_t const kND_nu = 0;
size_t const kND_293kA_nu = 0;
size_t const kND_280kA_nu = 1;
// The pairs of nonswap/numu and nueswap/nue are equivalent in terms of index
// position but not semantically identical. e.g. To make a full nue selection
// prediction you need the non swap (intrinsic) and nueswap (appeared files).
// However, both semantic meanings will not be used in the same list of objects
// and so can be safetly mapped onto the same indices.
size_t const kFD_nu_nonswap = 2;
size_t const kFD_nu_numu = 2;
size_t const kFD_nu_nueswap = 3;
size_t const kFD_nu_nue = 3;
size_t const kND_nub = 4;
size_t const kND_293kA_nub = 4;
size_t const kND_280kA_nub = 5;
size_t const kFD_nub_nonswap = 6;
size_t const kFD_nub_numu = 6;
size_t const kFD_nub_nueswap = 7;
size_t const kFD_nub_nue = 7;

size_t const kNPRISMConfigs_nu = 4;
size_t const kNPRISMConfigs = 8;
size_t const kNPRISMNDConfigs_nu = 2;
size_t const kNPRISMFDConfigs_nu = 2;
size_t const kNPRISMFDConfigs = 4;

inline bool IsNuConfig(size_t conf) { return conf < kNPRISMConfigs_nu; }

inline bool IsNDConfig(size_t conf) {
  return ((conf % kNPRISMConfigs_nu) <= 1);
}

inline bool IsND293kAConfig(size_t conf) {
  return ((conf % kNPRISMConfigs_nu) == 0);
}
inline bool IsND280kAConfig(size_t conf) {
  return ((conf % kNPRISMConfigs_nu) == 1);
}

inline size_t GetFDConfig(size_t conf) {
  if (IsNDConfig(conf)) {
    std::cout << "[ERROR]: Tried to get FD sub-config from an ND config("
              << conf << ")." << std::endl;
    abort();
  }
  return (((conf % kNPRISMConfigs_nu) - kNPRISMNDConfigs_nu) +
          (!IsNuConfig(conf) * kNPRISMFDConfigs_nu));
}

inline bool IsNuFDConfig(size_t conf) { return conf < kNPRISMFDConfigs_nu; }

inline size_t GetConfigFromFD(size_t conf) {
  if (conf >= kNPRISMFDConfigs) {
    std::cout << "[ERROR]: Tried to get config from an FD config(" << conf
              << ")." << std::endl;
    abort();
  }
  return (((conf % kNPRISMFDConfigs_nu) + kNPRISMNDConfigs_nu) +
          (!IsNuFDConfig(conf) * kNPRISMConfigs_nu));
}

inline bool IsNumuConfig(size_t conf) {
  return IsNDConfig(conf) || ((conf % kNPRISMConfigs_nu) == kFD_nu_numu);
}

inline bool IsNueConfig(size_t conf) {
  return !IsNDConfig(conf) && ((conf % kNPRISMConfigs_nu) == kFD_nu_nue);
}

inline size_t GetConfigNueSwap(size_t conf) {
  if (IsNDConfig(conf)) {
    std::cout << "[ERROR]: Tried to get FD sub-config from an ND config("
              << conf << ")." << std::endl;
    abort();
  }
  if (IsNueConfig(conf)) {
    return conf;
  }
  // Is a numu config
  return conf + 1;
}

inline size_t GetConfigNonSwap(size_t conf) {
  if (IsNDConfig(conf)) {
    std::cout << "[ERROR]: Tried to get FD sub-config from an ND config("
              << conf << ")." << std::endl;
    abort();
  }
  if (IsNueConfig(conf)) {
    return conf - 1;
  }
  // Is a numu config
  return conf;
}

inline std::string DescribeConfig(size_t conf) {
  switch (conf) {
  case kND_293kA_nu: {
    return "ND_293kA_nu";
  }
  case kND_280kA_nu: {
    return "ND_280kA_nu";
  }
  case kFD_nu_nonswap: {
    return "FD_nu_numu";
  }
  case kFD_nu_nueswap: {
    return "FD_nu_nue";
  }
  case kND_293kA_nub: {
    return "ND_293kA_nub";
  }
  case kND_280kA_nub: {
    return "ND_280kA_nub";
  }
  case kFD_nub_nonswap: {
    return "FD_nub_numu";
  }
  case kFD_nub_nueswap: {
    return "FD_nub_nue";
  }
  }
  std::cout << "Invalid conf: " << conf << std::endl;
  abort();
}

inline std::string DescribeFDConfig(size_t conf) {
  switch (GetConfigFromFD(conf)) {
  case kFD_nu_nonswap: {
    return "FD_nu_numu";
  }
  case kFD_nu_nueswap: {
    return "FD_nu_nue";
  }
  case kFD_nub_nonswap: {
    return "FD_nub_numu";
  }
  case kFD_nub_nueswap: {
    return "FD_nub_nue";
  }
  }
  std::cout << "Invalid FD conf: " << conf
            << " = (Determined FD conf: " << GetConfigFromFD(conf) << ")"
            << std::endl;
  abort();
}

inline size_t GetConfigFromNuChan(BeamChan nc, bool IsND) {
  if (IsND) {

    switch (nc.mode) {
    case BeamMode::kNuMode: {
      return kND_nu;
    }
    case BeamMode::kNuBarMode: {
      return kND_nub;
    }
    case BeamMode::kNuMode_293kA: {
      return kND_293kA_nu;
    }
    case BeamMode::kNuMode_280kA: {
      return kND_280kA_nu;
    }
    case BeamMode::kNuBarMode_293kA: {
      return kND_293kA_nub;
    }
    case BeamMode::kNuBarMode_280kA: {
      return kND_280kA_nub;
    }
    default: {
      std::cout << "[ERROR]: Explicit Horn Current Beam Mode required when "
                   "determining ND config."
                << std::endl;
      abort();
    }
    }
  } else {
    if ((nc.mode == BeamMode::kNuMode) && (nc.chan & NuChan::kNumu)) {
      return kFD_nu_numu;
    } else if ((nc.mode == BeamMode::kNuBarMode) &&
               (nc.chan & NuChan::kNumuBar)) {
      return kFD_nub_numu;
    } else if ((nc.mode == BeamMode::kNuMode) && (nc.chan & NuChan::kNue)) {
      return kFD_nu_nue;
    } else if ((nc.mode == BeamMode::kNuBarMode) &&
               (nc.chan & NuChan::kNueBar)) {
      return kFD_nub_nue;
    } else {
      std::cout << "[ERROR]: Invalid beam mode and channel: " << nc.mode << ", "
                << nc.chan << std::endl;
      abort();
    }
  }
}

inline size_t GetFDConfigFromNuChan(BeamChan nc) {
  return GetFDConfig(GetConfigFromNuChan(nc, false));
}

inline int FluxSpeciesPDG(NuChan fps) {
  if (fps & NuChan::kNumu) {
    return 14;
  } else if (fps & NuChan::kNue) {
    return 12;
  } else if (fps & NuChan::kNumuBar) {
    return -14;
  } else if (fps & NuChan::kNueBar) {
    return -12;
  }
  std::cout << "[ERROR]: Invalid NuChan: " << fps << std::endl;
  abort();
}

#define LOUDASSERT(a)                                                          \
  {                                                                            \
    if (!(a)) {                                                                \
      std::cout << "ASSERT FAILED: " << #a << std::endl;                       \
    }                                                                          \
    assert(a);                                                                 \
  }

inline void TestConfigDefinitions() {
  LOUDASSERT(IsNuConfig(kND_293kA_nu));
  LOUDASSERT(IsNDConfig(kND_293kA_nu));
  LOUDASSERT(IsND293kAConfig(kND_293kA_nu));
  LOUDASSERT(!IsND280kAConfig(kND_293kA_nu));
  LOUDASSERT(IsNumuConfig(kND_293kA_nu));
  LOUDASSERT(!IsNueConfig(kND_293kA_nu));

  LOUDASSERT(IsNuConfig(kND_280kA_nu));
  LOUDASSERT(IsNDConfig(kND_280kA_nu));
  LOUDASSERT(!IsND293kAConfig(kND_280kA_nu));
  LOUDASSERT(IsND280kAConfig(kND_280kA_nu));
  LOUDASSERT(IsNumuConfig(kND_280kA_nu));
  LOUDASSERT(!IsNueConfig(kND_280kA_nu));

  LOUDASSERT(!IsNuConfig(kND_293kA_nub));
  LOUDASSERT(IsNDConfig(kND_293kA_nub));
  LOUDASSERT(IsND293kAConfig(kND_293kA_nub));
  LOUDASSERT(!IsND280kAConfig(kND_293kA_nub));
  LOUDASSERT(IsNumuConfig(kND_293kA_nub));
  LOUDASSERT(!IsNueConfig(kND_293kA_nub));

  LOUDASSERT(!IsNuConfig(kND_280kA_nub));
  LOUDASSERT(IsNDConfig(kND_280kA_nub));
  LOUDASSERT(!IsND293kAConfig(kND_280kA_nub));
  LOUDASSERT(IsND280kAConfig(kND_280kA_nub));
  LOUDASSERT(IsNumuConfig(kND_280kA_nub));
  LOUDASSERT(!IsNueConfig(kND_280kA_nub));

  LOUDASSERT(IsNuConfig(kFD_nu_numu));
  LOUDASSERT(!IsNDConfig(kFD_nu_numu));
  LOUDASSERT(GetFDConfig(kFD_nu_numu) == 0);
  LOUDASSERT(IsNuFDConfig(GetFDConfig(kFD_nu_numu)));
  LOUDASSERT(GetConfigFromFD(GetFDConfig(kFD_nu_numu)) == kFD_nu_numu);
  LOUDASSERT(IsNumuConfig(kFD_nu_numu));
  LOUDASSERT(!IsNueConfig(kFD_nu_numu));
  LOUDASSERT(GetConfigNueSwap(kFD_nu_numu) == kFD_nu_nue);
  LOUDASSERT(GetConfigNonSwap(kFD_nu_numu) == kFD_nu_numu);

  LOUDASSERT(!IsNuConfig(kFD_nub_numu));
  LOUDASSERT(!IsNDConfig(kFD_nub_numu));
  LOUDASSERT(GetFDConfig(kFD_nub_numu) == 2);
  LOUDASSERT(!IsNuFDConfig(GetFDConfig(kFD_nub_numu)));
  LOUDASSERT(GetConfigFromFD(GetFDConfig(kFD_nub_numu)) == kFD_nub_numu);
  LOUDASSERT(IsNumuConfig(kFD_nub_numu));
  LOUDASSERT(!IsNueConfig(kFD_nub_numu));
  LOUDASSERT(GetConfigNueSwap(kFD_nub_numu) == kFD_nub_nue);
  LOUDASSERT(GetConfigNonSwap(kFD_nub_numu) == kFD_nub_numu);

  LOUDASSERT(IsNuConfig(kFD_nu_nue));
  LOUDASSERT(!IsNDConfig(kFD_nu_nue));
  LOUDASSERT(GetFDConfig(kFD_nu_nue) == 1);
  LOUDASSERT(IsNuFDConfig(GetFDConfig(kFD_nu_nue)));
  LOUDASSERT(GetConfigFromFD(GetFDConfig(kFD_nu_nue)) == kFD_nu_nue);
  LOUDASSERT(!IsNumuConfig(kFD_nu_nue));
  LOUDASSERT(IsNueConfig(kFD_nu_nue));
  LOUDASSERT(GetConfigNueSwap(kFD_nu_nue) == kFD_nu_nue);
  LOUDASSERT(GetConfigNonSwap(kFD_nu_nue) == kFD_nu_numu);

  LOUDASSERT(!IsNuConfig(kFD_nub_nue));
  LOUDASSERT(!IsNDConfig(kFD_nub_nue));
  LOUDASSERT(GetFDConfig(kFD_nub_nue) == 3);
  LOUDASSERT(!IsNuFDConfig(GetFDConfig(kFD_nub_nue)));
  LOUDASSERT(GetConfigFromFD(GetFDConfig(kFD_nub_nue)) == kFD_nub_nue);
  LOUDASSERT(!IsNumuConfig(kFD_nub_nue));
  LOUDASSERT(IsNueConfig(kFD_nub_nue));
  LOUDASSERT(GetConfigNueSwap(kFD_nub_nue) == kFD_nub_nue);
  LOUDASSERT(GetConfigNonSwap(kFD_nub_nue) == kFD_nub_numu);
}

struct PRISMAxisBlob {
  ana::HistAxis XProjection;
  ana::HistAxis OffAxisPosition;
  ana::HistAxis OffAxis280kAPosition;
};

ana::HistAxis GetEventRateMatchAxes(std::string const &binning = "uniform");
PRISMAxisBlob GetPRISMAxes(std::string const &varname,
                           std::string const &xbinning = "default",
                           std::string const &oabinning = "default");

ana::HistAxis TrueObservable(std::string const &obsvarname = "EProxy",
                             std::string const &binning = "uniform");

extern const ana::Cut kETrueLT10GeV;

extern const ana::Cut kETrue8GeV;

extern const ana::Cut kERecoProxy8GeV;

extern const ana::Cut kIsOutOfTheDesert;

// Use to weight by Exposure
extern const ana::Var kRunPlanWeight;

extern const ana::Cut kCut280kARun;
extern const ana::Cut kSel280kARun;
extern const ana::Var kSpecHCRunWeight;

ana::Cut GetNDSignalCut(bool UseOnAxisSelection = false, bool isNuMode = true);
ana::Cut GetFDSignalCut(bool UseOnAxisSelection = false, bool isNuMode = true,
                        bool isNuMu = true);

ana::Var GetNDWeight(std::string const &eweight = "", bool isNuMode = true);
ana::Var GetFDWeight(std::string const &eweight = "", bool isNuMode = true);

double Get280kAWeight_numu(double enu, bool isNu);
ana::Var GetNDSpecialRun(std::string const &SRDescriptor = "");

} // namespace PRISM

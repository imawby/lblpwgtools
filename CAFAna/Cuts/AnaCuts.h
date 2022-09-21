#pragma once

#include <cassert>
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/SRProxy.h"

namespace ana
{
    // same as the cvn!!
    inline bool IsInFV(const caf::SRProxy* sr)
    {
        const double minX(-360.0 + 50.0), maxX(360.0 - 50.0);
        const double minY(-600.0 + 50.0), maxY(600.0 - 50.0);
        const double minZ(50.0), maxZ(1394.0 - 150.0);

        if ((sr->recoVertex_x < minX) || (sr->recoVertex_x > maxX))
            return false;

        if ((sr->recoVertex_y < minY) || (sr->recoVertex_y > maxY))
            return false;

        if ((sr->recoVertex_z < minZ) || (sr->recoVertex_z > maxZ))
            return false;

        return true;
    }

    const Cut kPassFD_CVN_NUE([](const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        return (sr->cvnnue > 0.85 && sr->cvnnumu < 0.5);
    });

    const Cut kPassFD_CVN_NUMU([](const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        return (sr->cvnnumu > 0.5 && sr->cvnnue < 0.85);
    });

  const Cut kPassND_FHC_NUMU(
                  [](const caf::SRProxy* sr)
                  {
                    return (
			    sr->reco_numu && 
			    (sr->muon_contained || sr->muon_tracker) &&
			    sr->reco_q == -1 && 
			    sr->Ehad_veto<30);
		      });

    const Cut kPassND_RHC_NUMU(
                  [](const caf::SRProxy* sr)
                  {
                    return (
			    sr->reco_numu && 
			    (sr->muon_contained || sr->muon_tracker) &&
			    sr->reco_q == +1 && 
			    sr->Ehad_veto<30);
                  });

    inline bool PassNueSelectionFHC(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
	  return false;

        if (sr->selTrackPandizzleScore > sr->nuePandizzleCut)
            return false;

        if ((sr->selShowerPandrizzleScore < sr->nuePandrizzleCut) && (sr->selShowerJamPandrizzleScore < sr->nueJamPandrizzleCut))
            return false;

        return true;
    }

    inline bool PassNumuSelectionFHC(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        if (sr->selTrackPandizzleScore < sr->numuPandizzleCut)
            return false;

        return true;
    }

    inline bool PassNueSelectionRHC(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
	  return false;

        if (sr->selTrackPandizzleScore > sr->anuePandizzleCut)
            return false;

        if ((sr->selShowerPandrizzleScore < sr->anuePandrizzleCut) && (sr->selShowerJamPandrizzleScore < sr->anueJamPandrizzleCut))
            return false;

        return true;
    }

    inline bool PassNumuSelectionRHC(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        if (sr->selTrackPandizzleScore < sr->anumuPandizzleCut)
            return false;

        return true;
    }

    const Cut kIsNueSelectedFHC([](const caf::SRProxy* sr)
    {
        return PassNueSelectionFHC(sr);
    });

    const Cut kIsNumuSelectedFHC([](const caf::SRProxy* sr)
    {
        return (!PassNueSelectionFHC(sr) && PassNumuSelectionFHC(sr));
    });

    const Cut kIsNueSelectedRHC([](const caf::SRProxy* sr)
    {
        return PassNueSelectionRHC(sr);
    });

    const Cut kIsNumuSelectedRHC([](const caf::SRProxy* sr)
    {
        return (!PassNueSelectionRHC(sr) && PassNumuSelectionRHC(sr));
    });

}

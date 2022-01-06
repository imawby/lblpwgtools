#pragma once

#include <cassert>
#include "CAFAna/Core/Cut.h"
#include "StandardRecord/SRProxy.h"

namespace ana
{

  const Cut kPassFD_CVN_NUE(
                  [](const caf::SRProxy* sr)
                  {
                    return (sr->cvnnue > 0.85 && sr->cvnnumu < 0.5);
                  });

  const Cut kPassFD_CVN_NUMU(
                  [](const caf::SRProxy* sr)
                  {
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

    inline bool IsInFV(const caf::SRProxy* sr)
    {
        const double minX(-360.0 + 50.0), maxX(360.0 - 50.0);
        const double minY(-600.0 + 50.0), maxY(600.0 - 50.0);
        const double minZ(50.0), maxZ(1394.0 - 150.0);

        if ((sr->vtx_x < minX) || (sr->vtx_x > maxX))
            return false;

        if ((sr->vtx_y < minY) || (sr->vtx_y > maxY))
            return false;

        if ((sr->vtx_z < minZ) || (sr->vtx_z > maxZ))
            return false;

        return true;
    }

    inline bool PassNueSelection(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        if (sr->selTrackPandizzleScore > 0.72)
            return false;

        if (sr->selShowerPandrizzleScore < 0.42)
            return false;

        return true;
    }

    inline bool PassNumuSelection(const caf::SRProxy* sr)
    {
        if (!IsInFV(sr))
            return false;

        if (sr->selTrackPandizzleScore < 0.05)
            return false;

        return true;
    }

    const Cut kIsNueSelected([](const caf::SRProxy* sr)
    {
        return PassNueSelection(sr);
    });

    const Cut kIsNumuSelected([](const caf::SRProxy* sr)
    {
        return (!PassNueSelection(sr) && PassNumuSelection(sr));
    });

}

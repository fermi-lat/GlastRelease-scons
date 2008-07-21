#include "OnboardFilterTds/FilterStatus.h"

namespace OnboardFilterTds{

  float FilterStatus::getCalEnergy() const
  {
	  int energybits = m_stageEnergy & GFC_V2_STAGE_M_ENERGY;

	  // Drag the sign bit up and down
	  energybits <<= 8;
	  energybits >>= 8;

	  float energy = static_cast<float>(energybits) / 4.0;
    return energy;
  }


} // namespace

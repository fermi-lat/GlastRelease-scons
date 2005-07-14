#include "reconRootData/CalMipTrack.h"
#include <iostream>

ClassImp(CalMipTrack)

//-----------------------------------------------------------------------------------------------------------------
void CalMipTrack::print()
{
  std::cout << "---> writeOutPrint calMipTrack" << std::endl;
  std::cout << "---> point =" << getPoint().x() << " " << getPoint().y() << " " << getPoint().z() << std::endl;
  std::cout << "---> dir   =" << getDir().x()   << " " << getDir().y()   << " " << getDir().z()   << std::endl;
  std::cout << "---> d2C / d2Edge / calEdge / arcLen = " << m_d2C << " " << m_d2Edge << " " << m_calEdge << " "<< m_arcLen << std::endl;
  std::cout << "---> ecor / ecorRms / chi2 / erm = " << m_ecor << " " << m_ecorRms << " " << m_chi2 << " " << m_erm << std::endl;
 
}

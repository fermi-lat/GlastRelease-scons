#include "reconRootData/CalMipTrack.h"
#include <iostream>

ClassImp(CalMipTrack)

//-----------------------------------------------------------------------------------------------------------------
void CalMipTrack::print()
{
    std::cout<<"-----------------------------------------------------------------"<<std::endl;
    std::cout<<" Print in CalMipTrack s Class"<<std::endl;

    std::cout<< "length      =" << m_length       <<std::endl;
    std::cout<< "dTrack2C    ="<< m_dTrack2C      <<std::endl;
    std::cout<< "dTrack2Edge ="<< m_dTrack2Edge   <<std::endl;
    std::cout<< "energy      ="<< m_energy        <<std::endl;
    std::cout<< "Ndof        =" << m_ndofTrack    <<std::endl;
    std::cout<< "Ki2         =" << m_ki2Track     <<std::endl;
    std::cout<< "PosX        =" << getPoint().x() <<std::endl;
    std::cout<< "PosY        =" << getPoint().y() <<std::endl;
    std::cout<< "PosZ        =" << getPoint().z() <<std::endl;
    std::cout<< "DirX        =" << getDir().x()   <<std::endl;
    std::cout<< "DirY        =" << getDir().y()   <<std::endl;
    std::cout<< "DirZ        =" << getDir().z()   <<std::endl;

    int comptor=0;
}

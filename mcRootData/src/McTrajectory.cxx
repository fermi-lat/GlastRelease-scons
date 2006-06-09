// File and Version Information:
//      $Header$
//
// Description:
//  The McTrajectory class contains the step by step position information for
//  particles tracked through GLAST during the simulation stage. 
//                                                
// Author:  Tracy Usher                    

#include "mcRootData/McTrajectory.h"
#include "Riostream.h"
#include <commonRootData/RootDataUtil.h>

ClassImp(McTrajectory)

McTrajectory::McTrajectory() : m_mcParticle(0), m_charge(0)
{
    m_pointCol = new TObjArray();

    m_pointCol->Clear();
}

McTrajectory::~McTrajectory() 
{    
    m_pointCol->Delete();
    delete m_pointCol;
    m_pointCol = 0;

	Clear();
}

void McTrajectory::initialize(McParticle* mcPart, UInt_t charge) 
{
    m_mcParticle = mcPart;
    m_charge     = charge;

    return;
}


void McTrajectory::Clear(Option_t *) 
{
    m_mcParticle = 0;
    m_charge     = 0;

    if (m_pointCol) m_pointCol->Clear();

    return;
}

// dummy data, just for tests
void McTrajectory::Fake( Int_t /*ievent*/, UInt_t rank, Float_t randNum ) {

    Clear() ;
    
    Float_t f = Float_t(rank) ;
    Float_t fr = f*randNum ;

    setMcTrajectoryCharge((UInt_t)fr) ;
    addMcPoint(new TVector3(fr,fr,fr)) ;
    addMcPoint(new TVector3(fr*2.,fr*2.,fr*2.)) ;

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t McTrajectory::CompareInRange( const McTrajectory & ref, const std::string & name ) const {

    Bool_t result = true ;
    
    result = COMPARE_IN_RANGE(McTrajectoryCharge) && result ;



    const McParticle *myPart = getMcParticle();
    const McParticle *refPart = ref.getMcParticle();
    result = rootdatautil::CompareInRange(*myPart,*refPart,"McParticle") && result ;


    //result = COMPARE_IN_RANGE(McParticle) && result ;  need to pass reference

    result = rootdatautil::TObjArrayCompareInRange<TVector3>(getMcPointCol(),ref.getMcPointCol(),"McPoint") && result ;
      
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

void McTrajectory::Print(Option_t *) const 
{
/*
    using namespace std;
    TObject::Print(option);
    std::cout.precision(2);
    std::cout << "Run: " << m_runId << " Event: " << m_eventId << std::endl;
    std::cout << m_particleCol->GetEntries() << " McParticles" << std::endl;
    std::cout << m_positionHitCol->GetEntries() << " McPositionHits" << std::endl;
    std::cout << m_integratingHitCol->GetEntries() 
        << " McIntegratingHits" << std::endl;
*/
}


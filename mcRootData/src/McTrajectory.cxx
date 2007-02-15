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

#include "McObjectManager.h"

ClassImp(McTrajectory)

McTrajectory::McTrajectory() 
{
//    m_pointCol = new TObjArray();

    m_pointCol.Clear();
}

McTrajectory::~McTrajectory() 
{    
    m_pointCol.Clear();
    //delete m_pointCol;
    //m_pointCol = 0;

	Clear();
}

void McTrajectory::initialize() 
{
    return;
}


void* McTrajectory::operator new(size_t size)
{
    McTrajectory* temp = McObjectManager::getPointer()->getNewMcTrajectory();

    // Clear out the TobjArray
    temp->m_pointCol.Clear();
    //delete temp->m_pointCol;
    //temp->m_pointCol = 0;

    return temp;
}

void* McTrajectory::operator new(size_t size, void* vp)
{
    return vp;
}

McTrajectory& McTrajectory::operator =(const McTrajectory& rhs)
{
//    m_pointCol = rhs.m_pointCol;
    for(int idx = 0; idx < rhs.getMcPointCol()->GetEntries(); idx++)
        m_pointCol.Add(rhs.getMcPointCol()->At(idx));

    return *this;
}

void McTrajectory::Clear(Option_t *) 
{
//    if (m_pointCol) m_pointCol->Clear();
    m_pointCol.Clear();

    return;
}

// dummy data, just for tests
void McTrajectory::Fake( Int_t /*ievent*/, UInt_t rank, Float_t randNum ) {

    Clear() ;
    
    Float_t f = Float_t(rank) ;
    Float_t fr = f*randNum ;

    McTrajectoryPoint* point = new McTrajectoryPoint();
    point->Fake(0, rank, randNum);
    addMcPoint(point);

    point = new McTrajectoryPoint();
    point->Fake(0, 2*rank, randNum);
    addMcPoint(point);
}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t McTrajectory::CompareInRange( const McTrajectory & ref, const std::string & name ) const {

    Bool_t result = true ;
    
//    result = COMPARE_IN_RANGE(McTrajectoryCharge) && result 

//    const McParticle *myPart = getMcParticle();
//    const McParticle *refPart = ref.getMcParticle();
//    result = rootdatautil::CompareInRange(*myPart,*refPart,"McParticle") && result ;


    //result = COMPARE_IN_RANGE(McParticle) && result ;  need to pass reference

    result = rootdatautil::TObjArrayCompareInRange<McTrajectoryPoint>(getMcPointCol(),ref.getMcPointCol(),"McPoint") && result ;
      
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

ClassImp(McTrajectoryPoint)

McTrajectoryPoint::McTrajectoryPoint() 
{
    Clear();
}

void McTrajectoryPoint::initialize(VolumeIdentifier id, Float_t energy, const TVector3& position)
{
    m_volumeID = id;
    m_energy   = energy;
    m_point    = position;
}

void McTrajectoryPoint::Clear(Option_t *option)
{
    m_volumeID.Clear();
    m_energy = 0.;
    m_point = TVector3(0.,0.,0.);
}

void* McTrajectoryPoint::operator new(size_t size)
{
    McTrajectoryPoint* temp = McObjectManager::getPointer()->getNewMcTrajectoryPoint();

    return temp;
}

void* McTrajectoryPoint::operator new(size_t size, void* vp)
{
    return vp;
}

McTrajectoryPoint& McTrajectoryPoint::operator =(const McTrajectoryPoint& rhs)
{
    m_volumeID = rhs.m_volumeID;
    m_energy   = rhs.m_energy;
    m_point    = rhs.m_point;

    return *this;
}

void McTrajectoryPoint::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) 
{
    Clear() ;
    
    Float_t f = Float_t(rank) ;
    Float_t fr = f*randNum ;
            
    VolumeIdentifier id ;
    id.append(1) ;

    TVector3 position(fr,fr,fr);
    initialize(id,f,position);
}

Bool_t McTrajectoryPoint::CompareInRange( const McTrajectoryPoint & ref, const std::string & name ) const
{
    bool result = 0;

    result = COMPARE_IN_RANGE(VolumeID) && result ;
    result = COMPARE_IN_RANGE(Energy) && result ;
    result = COMPARE_IN_RANGE(Point) && result ;

    return result;
}

void McTrajectoryPoint::Print(Option_t *option) const
{
}

#include "reconRootData/TkrVecPointsLink.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrVecPointsLink)

TkrVecPointsLink::TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang)
{
    initializeInfo(firstPoint, secondPoint, 0, ang);
}

void TkrVecPointsLink::initializeInfo(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, UInt_t status, Double_t ang)
{
    m_firstVecPoint  = firstPoint;
    m_secondVecPoint = secondPoint;
    m_statusBits     = status;
    m_maxScatAngle   = ang;

    // Are start and end points in the same tower?
    if (firstPoint->getXCluster()->getTkrId().getTowerX() == secondPoint->getXCluster()->getTkrId().getTowerX() &&
        firstPoint->getXCluster()->getTkrId().getTowerY() == secondPoint->getXCluster()->getTkrId().getTowerY())
        m_statusBits |= SAMETOWER;

    // Get the slopes in both X and Y to make the direction of this link
    double deltaX  = firstPoint->getXCluster()->getPosition().x() - secondPoint->getXCluster()->getPosition().x();
    double deltaXZ = firstPoint->getXCluster()->getPosition().z() - secondPoint->getXCluster()->getPosition().z();
    double slopeX  = deltaX / deltaXZ;

    double deltaY  = firstPoint->getYCluster()->getPosition().y() - secondPoint->getYCluster()->getPosition().y();
    double deltaYZ = firstPoint->getYCluster()->getPosition().z() - secondPoint->getYCluster()->getPosition().z();
    double slopeY  = deltaY / deltaYZ;

    TVector3 init_dir(-slopeX, -slopeY, -1.);
    m_vector       = init_dir.Unit();

    // Now get the position which we'll set to be the bilayer z, between the two clusters
    double linkZ   = 0.5 * (firstPoint->getXCluster()->getPosition().z() + firstPoint->getYCluster()->getPosition().z());
    double linkX   = firstPoint->getXCluster()->getPosition().x() + slopeX * (linkZ - firstPoint->getXCluster()->getPosition().z());
    double linkY   = firstPoint->getYCluster()->getPosition().y() + slopeY * (linkZ - firstPoint->getYCluster()->getPosition().z());

    m_position     = TVector3(linkX, linkY, linkZ);
}
    
const TVector3 TkrVecPointsLink::getPosition(double z) const
{
    // Set the position to the links position
    TVector3 position = m_position;

    // First get the arclength of the step
    double arcLen = (z - position.z()) / m_vector.z();

    // Advance the link position to the new z
    position = position + arcLen * m_vector;

    return position;
}

void TkrVecPointsLink::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_firstVecPoint   = 0;
    m_secondVecPoint  = 0;
    m_statusBits      = 0;
    m_position        = TVector3(0.,0.,0.);
    m_vector          = TVector3(0.,0.,0.);
    m_maxScatAngle    = 0.;
}

void TkrVecPointsLink::Print(Option_t *option) const {
    TObject::Print(option);
}



//======================================================
// For Unit Tests
//======================================================

void TkrVecPointsLink::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrVecPointsLink::CompareInRange( const TkrVecPointsLink & /*ref*/, const std::string & name ) const {

    bool result = true ;

    //result = COMPARE_IN_RANGE(Energy) && result ;

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

void* TkrVecPointsLink::operator new(size_t /*size*/)
{
    TkrVecPointsLink* temp = ReconObjectManager::getPointer()->getNewTkrVecPointsLink();

    return temp;
}

void* TkrVecPointsLink::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrVecPointsLink::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}


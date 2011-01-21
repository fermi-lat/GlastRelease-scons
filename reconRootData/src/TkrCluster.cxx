#include "reconRootData/TkrCluster.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrCluster)

TkrCluster::TkrCluster() {
    Clear();
}

TkrCluster::TkrCluster(commonRootData::TkrId tkrId, UInt_t istrip0, UInt_t istripf, 
                       const TVector3& position, UInt_t rawToT, Float_t ToT, UInt_t status, UInt_t nBad)
{
        // Purpose and method: makes a cluster with attributes
        // Input:  tkrId is a valid TkrId which identifies where the cluster is
    //         strip0 is the first strip in the cluster
    //         stripf is the last strip in the cluster
    //         position is the position of the cluster
    //         ToT is the time-over-threshold for this cluster
        //         flag tells if the cluster is "in use" or not
    //         id is the sequential cluster number

        Clear();
        
    m_tkrId    = tkrId;
    m_strip0   = istrip0;
    m_stripf   = istripf;        
    m_nBad     = nBad;
    m_position = position;
    m_rawToT   = rawToT;
    m_ToT      = ToT;        
    m_status   = status;
}


TkrCluster::~TkrCluster() 
{
    Clear();
}



void TkrCluster::Clear(Option_t* /* option */) 
{
    m_tkrId    = commonRootData::TkrId();
    m_status   = 0;
    m_strip0   = 0;
    m_stripf   = 0;
    m_ToT      = 0.0;
    m_position = TVector3(0.0, 0.0, 0.0);
    m_nBad     = 0;
}

void TkrCluster::Print(Option_t *option) const 
{
    TObject::Print(option);
    using namespace std;

        cout << "(TkrId): (" << m_tkrId << endl;
        cout << "First Strip, Last Strip " << m_strip0 << ", " << m_stripf << endl;
        cout << "ToT, position(x, y, z) " << m_ToT  << " ("
                << m_position.X() << ", " << m_position.Y() << ", " << m_position.Z() 
                << ")" << endl;
}


//======================================================
// For Unit Tests
//======================================================

void TkrCluster::Fake( Int_t ievent, UInt_t rank, Float_t randNum ) {

    Float_t f = Float_t(ievent) ;
    
    m_tkrId.Fake(ievent,rank,randNum) ;  
    m_strip0 = 4 ;
    m_stripf = 10 ;
    m_nBad = 0 ;
    m_position = TVector3(f,2.*f,3.*f) ;
    m_rawToT = 50 ;
    m_ToT = f ;    
    m_status = 1 ;

}

#define COMPARE_IN_RANGE(att,text) rootdatautil::CompareInRange(att,ref.att,text)

Bool_t TkrCluster::CompareInRange( const TkrCluster & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(m_tkrId,"volume id") && result ;
    result = COMPARE_IN_RANGE(m_strip0,"Initial Strip Address") && result ;
    result = COMPARE_IN_RANGE(m_stripf,"Final Strip Address") && result ;
    result = COMPARE_IN_RANGE(m_nBad,"Number of Bad Strips")  && result ;
    result = COMPARE_IN_RANGE(m_position,"Space Position")  && result ;
    result = COMPARE_IN_RANGE(m_rawToT,"Raw ToT")  && result ;
    result = COMPARE_IN_RANGE(m_ToT,"ToT") && result ;
    result = COMPARE_IN_RANGE(m_status,"Status") && result ;
    result = COMPARE_IN_RANGE(getStrip(),"Center") && result ; // redundant
    result = COMPARE_IN_RANGE(getSize(),"Size") && result ; // redundant
    result = COMPARE_IN_RANGE(hitFlagged(),"Flag") && result ; // redundant

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

void* TkrCluster::operator new(size_t /*size*/)
{
    TkrCluster* temp = ReconObjectManager::getPointer()->getNewTkrCluster();

    // Since we recycle, make sure these member functions are cleared
    //temp->m_hitCol.clear();

    return temp;
}

void* TkrCluster::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrCluster::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}


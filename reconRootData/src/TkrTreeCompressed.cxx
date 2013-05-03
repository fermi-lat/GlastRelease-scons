#include "reconRootData/TkrTreeCompressed.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrTreeCompressed)

TkrTreeCompressed::TkrTreeCompressed() : m_headNode(0), 
                     m_bestLeaf(0), 
                     m_secondLeaf(0), 
                     m_axisParams(0), 
                     m_bestBranchAngleToAxis(0.), 
                     m_axisSeededAngleToAxis(0.)
{}

void TkrTreeCompressed::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_headNode              = 0;
    m_bestLeaf              = 0;
    m_secondLeaf            = 0;
    m_bestBranchAngleToAxis = 0.;
    m_axisSeededAngleToAxis = 0.;

    if (m_axisParams) delete m_axisParams;
    m_axisParams = 0;

    TObjArray::Clear();
}

void TkrTreeCompressed::Print(Option_t *option) const {
    TObject::Print(option);
}

void TkrTreeCompressed::initializeInfo(TkrVecNodeCompressed* head, 
                                       TkrVecNodeCompressed* best, 
                                       TkrVecNodeCompressed* second, 
                                       TkrFilterParams*      filter,
                                       Double_t              bestBranchAngToAxis,
                                       Double_t              axisSeededAngToAxis)
{
    m_headNode              = head;
    m_bestLeaf              = best;
    m_secondLeaf            = second;
    m_axisParams            = filter;
    m_bestBranchAngleToAxis = bestBranchAngToAxis;
    m_axisSeededAngleToAxis = axisSeededAngToAxis;

    return;
}


//======================================================
// For Unit Tests
//======================================================

void TkrTreeCompressed::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTreeCompressed::CompareInRange( const TkrTreeCompressed & /*ref*/, const std::string & name ) const {

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

void* TkrTreeCompressed::operator new(size_t /*size*/)
{
    TkrTreeCompressed* temp = ReconObjectManager::getPointer()->getNewTkrTreeCompressed();

    return temp;
}

void* TkrTreeCompressed::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrTreeCompressed::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}


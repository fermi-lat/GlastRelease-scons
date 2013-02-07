#include "reconRootData/TkrTree.h"
#include "Riostream.h"

#include "ReconObjectManager.h"

ClassImp(TkrTree)

TkrTree::TkrTree() : m_headNode(0), 
                     m_bestLeaf(0), 
                     m_secondLeaf(0), 
                     m_axisParams(0), 
                     m_bestBranchAngleToAxis(0.), 
                     m_axisSeededAngleToAxis(0.)
{}

void TkrTree::Clear(Option_t* /* option */) {
    // Make sure everything is properly zeroed
    m_headNode              = 0;
    m_bestLeaf              = 0;
    m_secondLeaf            = 0;
    m_bestBranchAngleToAxis = 0.;
    m_axisSeededAngleToAxis = 0.;

    if (m_axisParams) delete m_axisParams;
    m_axisParams = 0;

    TRefArray::Clear();
}

void TkrTree::Print(Option_t *option) const {
    TObject::Print(option);
}

void TkrTree::initializeInfo(TkrVecNode*      head, 
                             TkrVecNode*      best, 
                             TkrVecNode*      second, 
                             TkrFilterParams* filter,
                             Double_t         bestBranchAngToAxis,
                             Double_t         axisSeededAngToAxis)
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

void TkrTree::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTree::CompareInRange( const TkrTree & /*ref*/, const std::string & name ) const {

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

void* TkrTree::operator new(size_t /*size*/)
{
    TkrTree* temp = ReconObjectManager::getPointer()->getNewTkrTree();

    return temp;
}

void* TkrTree::operator new(size_t /*size*/, void* vp)
{
    return vp;
}

void TkrTree::operator delete(void* /*p*/)
{
    // Since we let ReconObjectManager handles memory, nothing to do here
    return;
}


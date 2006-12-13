#include "gcrSelectRootData/GcrSelect.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrSelect)

GcrSelect::GcrSelect() 
{
    m_gcrSelectedXtalCol = 0;
    m_gcrSelectVals = 0;
    m_indGcrSelectXtal =-1;
}

GcrSelect::~GcrSelect() 
{
    
    if (m_gcrSelectedXtalCol) 
    {
        m_gcrSelectedXtalCol->Delete();
        delete m_gcrSelectedXtalCol;
        m_gcrSelectedXtalCol = 0;
    }
    if(m_gcrSelectVals){m_gcrSelectVals = 0;}

}

void GcrSelect::initialize() 
{
    m_indGcrSelectXtal =-1;
    //if (!m_gcrSelectedXtalCol) m_gcrSelectedXtalCol = new TObjArray();

}

void GcrSelect::Clear(Option_t* /* option */) 
{
    if (m_gcrSelectedXtalCol) m_gcrSelectedXtalCol->Clear("C");
    m_indGcrSelectXtal = -1;
    if (m_gcrSelectVals) m_gcrSelectVals =0;
    
}

GcrSelectedXtal* GcrSelect::addGcrSelectedXtal() {

    if (!m_gcrSelectedXtalCol) m_gcrSelectedXtalCol = new TClonesArray("GcrSelectedXtal", 1);
    ++m_indGcrSelectXtal;
    TClonesArray &localCol = *m_gcrSelectedXtalCol;
    new(localCol[m_indGcrSelectXtal]) GcrSelectedXtal();
    return ((GcrSelectedXtal*) (localCol[m_indGcrSelectXtal]));

}

void GcrSelect::Print(Option_t *option) const 
{
    TObject::Print(option);
    std::cout
       << " # gcrSelectedXtal : " << m_indGcrSelectXtal+1
	<< std::endl;
}


//======================================================
// For Unit Tests
//======================================================

void GcrSelect::Fake( Int_t ievent, Float_t randNum ) {

    GcrSelectVals* m_gcrSelectVals = new GcrSelectVals();
    m_gcrSelectVals->Fake(ievent, randNum);
    addGcrSelectVals(m_gcrSelectVals);

    const Int_t numXtals = 5;
    Int_t ixtal;
    for (ixtal = 0; ixtal < numXtals; ixtal++) {
        GcrSelectedXtal* xtal = addGcrSelectedXtal();
        xtal->Fake(ievent, ixtal, randNum);

    }
}

#define COMPARE_TOBJ_ARRAY_IN_RANGE(T,m) rootdatautil::TObjArrayCompareInRange<T>(m,ref.m)

Bool_t GcrSelect::CompareInRange(const GcrSelect & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_TOBJ_ARRAY_IN_RANGE(GcrSelectedXtal,getGcrSelectedXtalCol()) && result ;

    result = ((GcrSelectVals*)m_gcrSelectVals)->CompareInRange(*(GcrSelectVals*)(ref.getGcrSelectVals()), "GcrSelectVals") && result;

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

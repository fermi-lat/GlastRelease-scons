#include "reconRootData/CalRecon.h"
#include "gcrSelectRootData/GcrSelect.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrSelect)

GcrSelect::GcrSelect() 
{
    //m_gcrXtalCol = 0;
    m_gcrSelectedXtalCol = 0;
    m_gcrSelectVals = 0;
}

GcrSelect::~GcrSelect() 
{
    
    /**if (m_gcrXtalCol) 
    {
        m_gcrXtalCol->Delete();
        delete m_gcrXtalCol;
        m_gcrXtalCol = 0;
    }*/
    
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
    //if (!m_gcrXtalCol) m_gcrXtalCol = new TObjArray();
    if (!m_gcrSelectedXtalCol) m_gcrSelectedXtalCol = new TObjArray();
    //if (!m_gcrSelectVals) m_gcrSelectVals = new TObject();

}

void GcrSelect::Clear(Option_t* /* option */) 
{
    //if (m_gcrXtalCol) m_gcrXtalCol->Delete();
    if (m_gcrSelectedXtalCol) m_gcrSelectedXtalCol->Delete();
    if (m_gcrSelectVals) m_gcrSelectVals =0;
    
}

void GcrSelect::Print(Option_t *option) const 
{
    TObject::Print(option);
    std::cout
        //<< " # gcrXtal : " << m_gcrXtalCol->GetEntries() 
       << " # gcrSelectedXtal : " << m_gcrSelectedXtalCol->GetEntries() 
	
	
	<< std::endl;
}


//======================================================
// For Unit Tests
//======================================================

/**void CalRecon::Fake( Int_t ievent, Float_t randNum ) {

    const UInt_t NUM_XTALS = 10;
    const UInt_t NUM_EVENT_ENERGIES = 2 ;
    const UInt_t NUM_CLUSTERS = 20;
    const UInt_t NUM_MIP_TRACKS = 5;

    initialize() ;
    Clear() ;

    UInt_t ienergy;
    for (ienergy = 0; ienergy < NUM_EVENT_ENERGIES ; ienergy++ ) {
        CalEventEnergy * energy = new CalEventEnergy() ;
        energy->Fake(ievent,ienergy,randNum) ;
        addCalEventEnergy(energy);
    }

    UInt_t icluster;
    for (icluster = 0; icluster < NUM_CLUSTERS ; icluster++ ) {
        CalCluster *cluster = new CalCluster();
        cluster->Fake(ievent,icluster,randNum) ;
        addCalCluster(cluster);
    }

    UInt_t ixtal;
    for (ixtal = 0; ixtal < NUM_XTALS ; ixtal++) {
        CalXtalRecData *xtal = new CalXtalRecData();
        xtal->Fake(ievent,ixtal,randNum) ;
        addXtalRecData(xtal);
    }

    unsigned int imip;
    for (imip = 0; imip < NUM_MIP_TRACKS ; imip++) {
        CalMipTrack * mipTrack = new CalMipTrack ;
        mipTrack->Fake(ievent,imip,randNum) ;
        addCalMipTrack(mipTrack);

    }

}

#define COMPARE_TOBJ_ARRAY_IN_RANGE(T,m) rootdatautil::TObjArrayCompareInRange<T>(m,ref.m)

Bool_t CalRecon::CompareInRange( CalRecon & ref, const std::string & name ) {

    bool result = true ;

    result = COMPARE_TOBJ_ARRAY_IN_RANGE(CalEventEnergy,getCalEventEnergyCol()) && result ;
    result = COMPARE_TOBJ_ARRAY_IN_RANGE(CalCluster,getCalClusterCol()) && result ;
    result = COMPARE_TOBJ_ARRAY_IN_RANGE(CalXtalRecData,getCalXtalRecCol()) && result ;
    result = COMPARE_TOBJ_ARRAY_IN_RANGE(CalMipTrack,getCalMipTrackCol()) && result ;

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
*/
//C.L. for DEBUG, overload of addGcrXtal method:
/**void GcrSelect::addGcrSelectedXtal(GcrSelectedXtal* gcrSelectedXtal)     { 
    
      
      std::cout<<"GcrSelect::addGcrSelectedXtal BEGIN"<<std::endl ;
      std::cout<< "gcrSelectedXtal->getXtalId()= " << gcrSelectedXtal->getXtalId().getTower() << "/" << gcrSelectedXtal->getXtalId().getLayer() << "/" << gcrSelectedXtal->getXtalId().getColumn() 
          << "gcrSelectedXtal-> getRowEnergy=" << gcrSelectedXtal-> getRowEnergy() 
          << "gcrSelectedXtal-> getCorrEnergy=" << gcrSelectedXtal-> getCorrEnergy() 
          << "gcrSelectedXtal-> getPathLength=" << gcrSelectedXtal-> getPathLength() 
          << "gcrSelectedXtal-> getInferedZ=" << gcrSelectedXtal-> getInferedZ() 
	  << "gcrSelectedXtal-> getCrossedFaces=" << gcrSelectedXtal-> getCrossedFaces()
	  << "gcrSelectedXtal-> getEntryPoint.x=" << gcrSelectedXtal-> getEntryPoint().x()
	  
	  << std::endl;   


      m_gcrSelectedXtalCol->Add(gcrSelectedXtal);
      
            std::cout<<"GcrSelect::addGcrSelectedXtal END"<<std::endl ;

    
    
     }*/


// $Header$
// Original author: Ian Gable
// GaudiGlastTuple.cxx: implementation of the GaudiGlastTuple class.
//
//////////////////////////////////////////////////////////////////////

#include "CalRecon/GaudiGlastTuple.h"

static std::string tupleRoot("/NTUPLES/");

/*! Here we need to pass in a pointer to the NTupleSvc, which you get from
    from the ntupleSvc() call back in the algorithm.
*/
GaudiGlastTuple::GaudiGlastTuple (INTupleSvc* ntSvc, const char* title, const char * logicalFileName ) 
:  m_ntSvc(ntSvc)
{
    m_fileName=tupleRoot+std::string(logicalFileName);

    NTupleDirPtr dir(m_ntSvc, m_fileName);
    if( dir ) {
        m_nt = new NTuplePtr( dir, "/1"); 
        if(!(*m_nt)){
            (*m_nt) = m_ntSvc->book(m_fileName,1,CLID_RowWiseTuple, title);
        }
    }
}

GaudiGlastTuple::~GaudiGlastTuple ()
{
    /// delete the NTuple::Item<float> objects;
    for(iterator tupit = m_ntupleItemList.begin();
                tupit != m_ntupleItemList.end(); ++tupit ) {
        delete (*tupit).first;
    }
    delete m_nt;
}

/*! 
*/
void GaudiGlastTuple::addItem (const char* name, const float *datum)
{
    // create an item for the Gaudi tuple, save pointer along with pointer to local data
    NTuple::Item<float>* ntItem = new NTuple::Item<float>;
    m_ntupleItemList.push_back( std::make_pair(ntItem,datum) );

    if(*m_nt) (*m_nt)->addItem( name ,*ntItem);
}


/*! This is where the information in the NTuple actually gets added
    at the end of the reconstruction.
*/
void GaudiGlastTuple::fill()
{
    for( iterator tupit = m_ntupleItemList.begin();
            tupit !=  m_ntupleItemList.end(); ++tupit) {
        *((*tupit).first) = *((*tupit).second);
    }
    StatusCode status = m_ntSvc->writeRecord(m_fileName+"/1");
}

/*! Need to dereference the pointer before you hand it back.
*/
NTuplePtr GaudiGlastTuple::getNTuple(){
    return *m_nt;
}


// GaudiGlastTuple.cxx: implementation of the GaudiGlastTuple class.
//
//////////////////////////////////////////////////////////////////////

#include "CalRecon\GaudiGlastTuple.h"

/*! Here we need to pass in a pointer to the NTupleSvc, which you get from
    from the ntupleSvc() call back in the algorithm.
*/
GaudiGlastTuple::GaudiGlastTuple (const char* title, INTupleSvc* ntSvc) :Tuple(title)
{


    m_ntSvc = ntSvc;
    m_ntName = title;

    m_file1 = new NTupleFilePtr(m_ntSvc, "/NTUPLES/FILE1");
    if((*m_file1)){
        m_nt = new NTuplePtr(m_ntSvc, "/NTUPLES/FILE1/1");
        if(!(*m_nt)){
            (*m_nt) = m_ntSvc->book("/NTUPLES/FILE1",1,CLID_ColumnWiseTuple, m_ntName);
        }
    }

    
}

/*! Needs to be filled in.
*/
GaudiGlastTuple::~GaudiGlastTuple ()
{
}

/*! Adds items to the m_float_array as well as directly to the.
    Ntuple
*/
void GaudiGlastTuple::addItem (const char* name, const float *datum)
{
    StatusCode status;
    m_ntupleItemList.push_back(new NTuple::Item<float>);
    m_float_array.push_back(datum);

    if(*m_nt){
       status = (*m_nt)->addItem( name ,*m_ntupleItemList[m_ntupleItemList.size()-1]);
       *m_ntupleItemList[m_ntupleItemList.size()-1] = *datum;
    }
    
}

void GaudiGlastTuple::dumpData ()
{
   //don't need anymore  m_column = begin();  // reset column pointer to begin
}

/*! This is where the information in the NTuple actually gets added
    at the end of the reconstruction.
*/
void GaudiGlastTuple::fill()
{
    StatusCode status;
    int i = 0;
    for(std::vector<const float*>::iterator it = m_float_array.begin();
                                            it != m_float_array.end();  
                                            ++it)
    {
	*m_ntupleItemList[i++] = **it;
    }
    
    // This is not working here for reasons unknown.
    //status = m_ntSvc->writeRecord("/NTUPLES/FILE1/1");
    
}

/*! Need to dereference the pointer before you hand it back.
*/
NTuplePtr GaudiGlastTuple::getNTuple(){
    return *m_nt;
}


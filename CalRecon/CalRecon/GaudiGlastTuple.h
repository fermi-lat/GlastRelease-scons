// $Header$
// Original author: Ian Gable
// GaudiGlastTuple.h: interface for the GaudiGlastTuple class.
//
//////////////////////////////////////////////////////////////////////

#ifndef GaudiGlastTuple_H
#define GaudiGlastTuple_H

#include <vector>

#include "Gaudi/Interfaces/INTupleSvc.h"
#include "Gaudi/Interfaces/INTuple.h"
#include "Gaudi/NTupleSvc/NTuple.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "Gaudi/Kernel/StatusCode.h"


/*! \class GaudiGlastTuple
\brief Wrapper for the Ntuple Servise so that this can plug into the SummaryData
template. Will create a row-wise ntuple.

*/
class GaudiGlastTuple {
    
public:
    /** 
       @param ntSvc The ntuple service
       @param title Title of the tuple
       @param logicalFileName For example, "FILE1". This must be defined by the job options, as
       NTupleSvc.Output = {"FILE1 DATAFILE='/NTUPLES/CalRecon.root' OPT='NEW'"};

  */
    GaudiGlastTuple(INTupleSvc* ntSvc, const char* title, const char* logicalFileName );
    ~GaudiGlastTuple();
    
    //! setup: add a new item, by name and a pointer to the actual data
    void addItem (const char* name, const float *datum);
    
    /// Must be called to write out a row.
    void fill();
    
    /// access to the NTuplePtr
    NTuplePtr getNTuple();
    
private:
    
    INTupleSvc* m_ntSvc;   /// the tuple service
    NTuplePtr* m_nt; /// pointer to the actual ntuple

    std::string m_fileName; /// the logical file name, e.g. "/NTUPLES/FILE1"

    ///maintain a list of pairs of pointers to the Gaudi Ntuple::Item, and our float
    typedef std::vector< std::pair<NTuple::Item<float>*,const float*> > container;
    typedef container::iterator iterator;
    container m_ntupleItemList; 
};


#endif // GaudiGlastTuple_H

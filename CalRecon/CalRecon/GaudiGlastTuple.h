// GaudiGlastTuple.h: interface for the GaudiGlastTuple class.
//
//////////////////////////////////////////////////////////////////////

#ifndef GaudiGlastTuple_H
#define GaudiGlastTuple_H

#include <vector>

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"
#include "reconstruction/analysis/Tuple.h"


/*! \class GaudiGlastTuple
\brief Wrapper for the Ntuple Servise so that this can plug into the SummaryData
        template.

  */
class GaudiGlastTuple : public Tuple {

public:
  //! Constructor of this form must be provided
    GaudiGlastTuple(const char* name, INTupleSvc* ntSvc);
    ~GaudiGlastTuple();

    void addItem (const char* name, const float *datum);
    void fill();
    void dumpData ();
    NTuplePtr getNTuple();
;
private:

  std::vector<NTuple::Item<float>*> m_ntupleItemList;
  INTupleSvc* m_ntSvc;
  const char* m_ntName;
  NTupleFilePtr* m_file1;
  NTuplePtr* m_nt;
  std::vector<const float*> m_float_array;
};


#endif // CalRecoAlg_H

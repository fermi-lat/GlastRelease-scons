// McVertex.h: interface for the McVertex class.
//
//////////////////////////////////////////////////////////////////////

#ifndef McVertexCnv_H 
#define McVertexCnv_H 1

// Include files
//#include "Gaudi/Kernel/Converter.h"
#include "src/EventSelector/BaseCnv.h"

#include "GlastEvent/TopLevel/MCEvent.h"


// GlastEvent for creating the McEvent stuff
#include "GlastEvent/MonteCarlo/McParticle.h"
#include "GlastEvent/MonteCarlo/McVertex.h"



extern const CLID& CLID_McVertexCnv;

// Abstract factory to create the converter
template <class TYPE> class CnvFactory;



// Original Author: T Burnett 

/*! Handles the conversion of McVertex currently uses the IRF
    information.
*/

class McVertexCnv : public BaseCnv   {
  friend class CnvFactory<BaseCnv>;

public:
  /// Standard Constructor
  McVertexCnv( ISvcLocator* svc );
  /// Standard Destructor
  virtual ~McVertexCnv();
  /// class ID
  static const CLID& classID();
  /// Update complete object
  virtual StatusCode updateObj(IOpaqueAddress* pAddress, DataObject* pObject);
  /// Create object
  virtual StatusCode createObj(IOpaqueAddress* pAddress, DataObject*& refpObject);

private:

    McVertex*       m_root; // 

    McVertexCol*   m_vlist;
    McParticleCol* m_plist;


};
#endif // McVertex_H

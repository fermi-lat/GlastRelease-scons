// $Header$
#ifndef RootTkrScaleCnv_h
#define RootTkrScaleCnv_h

/** @class RootTkrScaleCnv

  Converter from Root to TCDS Tkr (muon) charge scale class

  @author J. Bogart
*/
#include "RootTkrBaseCnv.h"
#include "CalibData/Tkr/TkrScale.h"

template <class TYPE> class CnvFactory;

class RootTkrScaleCnv : public RootTkrBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<RootTkrScaleCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();


protected:


  RootTkrScaleCnv(ISvcLocator* svcs);

  virtual ~RootTkrScaleCnv() {}       // most likely nothing to do 

  /**
     The only substantive function in this class.  Given a
     ROOT file containing muon charge scale information, make 
     the corresponding TDS object.
   */
  virtual StatusCode i_createObj (const std::string& fname,
                                  DataObject*& refpObject);
  /*  virtual StatusCode i_createObj(TObject* pRootObj,       
      DataObject*& refpObject);  */
private:
  StatusCode readUnis(TTree* tree, int iTow, CalibData::TkrScaleCol* col);

};
#endif

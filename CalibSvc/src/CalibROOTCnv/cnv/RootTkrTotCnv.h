// $Header$
#ifndef RootTkrTotCnv_h
#define RootTkrTotCnv_h

/** @class RootTkrTotCnv

  Converter from Root to TCDS Tkr (charge injection) ToT class

  @author J. Bogart
*/
#include "RootTkrBaseCnv.h"
#include "CalibData/Tkr/TkrTot.h"

template <class TYPE> class CnvFactory;

class RootTkrTotCnv : public RootTkrBaseCnv {

  /// Friend needed for instantiation
  friend class CnvFactory<RootTkrTotCnv>;
public:
  const CLID& objType() const;
  static const CLID& classID();


protected:


  RootTkrTotCnv(ISvcLocator* svcs);

  virtual ~RootTkrTotCnv() {}       // most likely nothing to do 

  /**
     The only substantive function in this class.  Given a
     ROOT file containing charge injection ToT information, make 
     the corresponding TDS object.
   */
  virtual StatusCode i_createObj (const std::string& fname,
                                  DataObject*& refpObject);
  /*  virtual StatusCode i_createObj(TObject* pRootObj,       
      DataObject*& refpObject);  */
private:
  StatusCode readUnis(TTree* tree, int iTow, CalibData::TkrTotCol* col);

};
#endif

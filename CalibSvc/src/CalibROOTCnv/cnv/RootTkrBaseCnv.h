// $Header$
#ifndef RootTkrBaseCnv_h
#define RootTkrBaseCnv_h

/** @class RootTkrBaseCnv 

  Base class for TKR calibration converters from ROOT files to TCDS.
  All such converters need to do certain things, which are
  handled here.  Methods common to *all* calibrations are in the
  base class RootBaseCnv

  @author J. Bogart
*/

#include "RootBaseCnv.h"
#include "CalibData/Tkr/TkrBase.h"
#include "idents/TkrId.h"                    // not strictly necessary
#include "commonRootData/idents/TkrId.h"

class TTree;

class RootTkrBaseCnv : public RootBaseCnv {
public:

  RootTkrBaseCnv(ISvcLocator* svc, const CLID& clid);

  virtual ~RootTkrBaseCnv() {};

protected:
  /**
       All Tkr ROOT files consist of a tree per tower, named
          Tower00, Tower01, ...Tower33
          ...except they don't all need to be present
       Find the tree, if any, corresponding to specified bay#.
       Assumes ROOT file has already been opened
  */
  TTree* findTower(unsigned bay);
  TTree* findTower(unsigned row, unsigned col);
  /**
     Read in generic part (tower object) of this tree and store in
     CalibData::TkrBase::TkrTower object
  */
  StatusCode readTower(TTree* tree, unsigned bay, CalibData::TkrBase* pCol);


  /**
     Convert commonRootData::TkrId to idents::TkrId.  towerX, towerY,  
     tray and botTop fields must be present and these are the only ones
     used.  id arg. should point to already-allocated idents::TkrId object
  */
  StatusCode convertId(const commonRootData::TkrId& rootId, 
                       idents::TkrId* id);

  bool checkTower(const idents::TkrId& id, int iTow);

  /**
     Given a pointer to a TDS object which can be cast to "our" type, fill
     in corresponding information in the corresponding root class

     @param pTDSObj   Pointer to tds object to be converted
     @param pRootObj  Pointer to destination root object

     no implementation for now for TKR

  */
  virtual StatusCode fillRoot(CalibData::CalibBase* /* pTDSObj */, 
                              TObject* /* pRootObj */) {
    return StatusCode::FAILURE; 
}

  /// For writing ..no implementation for now for TKR
  virtual StatusCode createRoot(const std::string& /* fname */, 
                                CalibData::CalibBase* /* pTDSObj */) {
    return StatusCode::FAILURE; 
  }


  // Invoke factory method to make a uni.  Derived classes aren't friends
  // of CalibData::TkrBase, so have to do it for them
  CalibData::UniBase* makeUni(CalibData::TkrBase* col) {
    return col->m_factory->makeUni();
  }
private:
  // Needs to be private since only RootTkrBaseCnv is a friend of TkrBase,
  // not derived classes
  //    Can't think of any private data we might need

};

#endif

// $Header$
#include "GaudiKernel/MsgStream.h"
#include "RootTkrBaseCnv.h"

#include "CalibData/Tkr/TkrBase.h"
#include "CalibSvc/IInstrumentName.h"
#include "calibRootData/Tkr/TkrTower.h"

#include "facilities/Util.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TIterator.h"
#include "TKey.h"
#include "TString.h"


RootTkrBaseCnv::RootTkrBaseCnv(ISvcLocator* svc, const CLID& clid) :
  RootBaseCnv(svc, clid) {}

// A place to put some Tkr-specific utilities

TTree* RootTkrBaseCnv::findTower(unsigned bay) {
  unsigned row = bay / TKRBASE_MAXCOL;
  unsigned col = bay % TKRBASE_MAXCOL;
  return findTower(row, col);
}

TTree* RootTkrBaseCnv::findTower(unsigned row, unsigned col) {
  if( (row > TKRBASE_MAXROW) || (col > TKRBASE_MAXCOL)) { // out of range
    return 0;    // throw exception?
  }
  if (!m_inFile) return 0;   // need an open file to search..throw exception?

  // Form tree name
  std::string treeString("Tower");
  std::string temp("");
  int iRow = row;
  facilities::Util::itoa(iRow, temp);
  treeString += temp;
  temp="";
  int iCol = col;
  facilities::Util::itoa(iCol, temp);
  treeString += temp;
  const char* treeName = treeString.c_str();

  // Search for key
  TIter nextTopLevelKey(m_inFile->GetListOfKeys());
  TKey *keyTopLevel;
  TTree* t=0; 

  while  ( (keyTopLevel=(TKey*)nextTopLevelKey()) ) {

    TString name(keyTopLevel->GetName());
    TString className(keyTopLevel->GetClassName());

    if ((name.CompareTo(treeName)==0) && (className.CompareTo("TTree")==0))  {
      // Found It
      t = (TTree*)m_inFile->Get(keyTopLevel->GetName());
      return t;
    }
  }
  return 0;                // not found
}
StatusCode RootTkrBaseCnv::readTower(TTree* pTree, unsigned iTow,
                                     CalibData::TkrBase* pCol) {
  CalibData::TkrBase::TkrTower* tdsTow = pCol->makeTower(iTow);

  //  calibRootData::TkrTower* rootTow = 0;
  calibRootData::TkrTower* rootTow = new calibRootData::TkrTower;
  TObject* pTObj = rootTow;
  StatusCode ret = readRootObj(pTree, 
                               std::string("calibRootData::TkrTower"), 
                               pTObj);
  if (ret != StatusCode::SUCCESS) return ret;

  tdsTow->m_iRow = rootTow->getRow();
  tdsTow->m_iCol = rootTow->getCol();    
  tdsTow->m_hwserial = std::string((const char *) (rootTow->getSerial()) );

  delete rootTow;    // is this necessary?

  return StatusCode::SUCCESS;
}



bool RootTkrBaseCnv::checkTower(const idents::TkrId& id, int iTow) {
  return (((unsigned) iTow/TKRBASE_MAXCOL == id.getTowerY())  &&
          ((unsigned) iTow%TKRBASE_MAXCOL == id.getTowerX())   );
}

// For now only handle towerX, towerY, tray and botTop fields. These
// are the only ones needed for TKR calibrations currently, probably
// forever.
StatusCode RootTkrBaseCnv::convertId(const commonRootData::TkrId& rootId, 
                                     idents::TkrId* id) {
  unsigned towerX;
  unsigned towerY;
  unsigned tray;
  bool botTop;

  if (!id) return StatusCode::FAILURE;

  if ((!rootId.hasTowerX() ) || (!rootId.hasTowerY()) ||
      (!rootId.hasTray() ) || (!rootId.hasBotTop())   ) {
    return StatusCode::FAILURE;
  }

  towerX = rootId.getTowerX();
  towerY = rootId.getTowerY();
  tray = rootId.getTray();
  botTop = (rootId.getBotTop() != 0);

  *id = idents::TkrId(towerX, towerY, tray, botTop);
  return StatusCode::SUCCESS;
}



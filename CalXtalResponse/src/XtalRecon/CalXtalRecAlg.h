#ifndef CalXtalRecAlg_h
#define CalXtalRecAlg_h 

// LOCAL INCLUDES
#include "CalXtalResponse/IXtalRecTool.h"

// GLAST INCLUDES
#include "Event/Digi/CalDigi.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "ntupleWriterSvc/INTupleWriterSvc.h"


// EXTLIB INCLUDES
#include "GaudiKernel/Algorithm.h"

// STD INCLUDES

/** @class CalXtalRecAlg
    @brief  Calorimeter crystal reconstruction algorithm

    This algorithm reconstructs energy and position in each calorimeter crystal.
    See CalXtalResponse/XtalRecTool

    @author           A.Chekhtman
    @author           Zach Fetwrell

    $Header$
*/
class CalXtalRecAlg : public Algorithm
{
 public:
  CalXtalRecAlg(const std::string& name, ISvcLocator* pSvcLocator);
  /// initialize internal data members.
  StatusCode initialize();
  ///  Reconstruct ene & pos for all xtal hits in event
  StatusCode execute();
  /// required by Gaudi Algorithm class
  StatusCode finalize();

 private:
  ///  function for setting pointers to the input and output data in Gaudi TDS
  StatusCode retrieve();
    
  /// pointer to input data collection in TDS
  Event::CalDigiCol* m_calDigiCol;

  /// pointer to the output data collection in TDS
  Event::CalXtalRecCol* m_calXtalRecCol;
    
  /// pointer to CalResponse tool for converting xtal digi info -> energy 
  IXtalRecTool *m_xtalRecTool;

  /// pointer to event Header (evtId, runId, etc...)
  Event::EventHeader* m_evtHdr;

  /// name of IXtalRecTool instantiation
  StringProperty m_recToolName;

  /// name of CalXtalRecTuple tree.  Default = "" (no file).
  StringProperty m_tupleName;
  /// optional name of CalXtalRecTuple file.  Defaul = "" (use ntupleWriterSvc default)
  StringProperty m_tupleFilename;
  /// store current entry for CalTuple
  CalTupleEntry m_tupleEntry;
  /// pointer to tupleWriterSvc
  INTupleWriterSvc *m_tupleWriterSvc;
};

#endif

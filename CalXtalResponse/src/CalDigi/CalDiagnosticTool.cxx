// $Header$

// Include files

/** @file
    @author Z.Fewtrell
*/


// LOCAL
#include "CalDiagnosticTool.h"
#include "CalXtalResponse/ICalSignalTool.h"
#include "CalXtalResponse/ICalTrigTool.h"
#include "../CalCalib/IPrecalcCalibTool.h"

// GLAST
#include "CalUtil/CalDiagnosticWord.h"

// EXTLIB
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

// STD

using namespace CalUtil;
using namespace std;


//static ToolFactory<CalDiagnosticTool> s_factory;
//const IToolFactory& CalDiagnosticToolFactory = s_factory;
DECLARE_TOOL_FACTORY(CalDiagnosticTool);


CalDiagnosticTool::CalDiagnosticTool( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent)
  : AlgTool(type,name,parent),
    m_precalcCalibTool(0),
    m_calSignalTool(0),
    m_calTrigTool(0)
{
  declareInterface<ICalDiagnosticTool>(this);

  declareProperty("PrecalcCalibTool", m_precalcCalibName = "PrecalcCalibTool");
  declareProperty("CalSignalToolName", m_calSignalToolName = "CalSignalTool");
  declareProperty("CalTrigToolName", m_calTrigToolName = "CalTrigTool");

}

StatusCode CalDiagnosticTool::initialize() {
  MsgStream msglog(msgSvc(), name());   
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc;

  //-- jobOptions --//
  if ((sc = setProperties()).isFailure()) {
    msglog << MSG::ERROR << "Failed to set properties" << endreq;
    return sc;
  }

  // this tool may also be shared by other tools, global ownership
  sc = toolSvc()->retrieveTool("PrecalcCalibTool", 
                               m_precalcCalibName, 
                               m_precalcCalibTool,
                               0); // shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << m_precalcCalibName << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalSignalTool",
                               m_calSignalToolName,
                               m_calSignalTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << m_calSignalToolName << endreq;
    return sc;
  }

  sc = toolSvc()->retrieveTool("CalTrigTool",
                               m_calTrigToolName,
                               m_calTrigTool,
                               0); // intended to be shared
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  can't create " << m_calTrigToolName << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}


auto_ptr<LdfEvent::CalDiagnosticData> CalDiagnosticTool::getDiagnosticData(const CalUtil::TwrNum twr,
                                                                           const CalUtil::LyrNum lyr) {
  /// initialize return value to null pointer.
  auto_ptr<LdfEvent::CalDiagnosticData> calDiagData;

  /// cal diag contains 2 trigger bits (one per diode size)
  CalDiagnosticWord::CalDiagTrigBits trigBits;
  if (calcTrigBits(twr, lyr, trigBits).isFailure())
    return calDiagData; ///< return NULL ptr

  /// cal diag contains 12 lac bits (one per xtal column)
  CalDiagnosticWord::CalDiagLACBits lacBits;
  if (calcLACBits(twr, lyr, lacBits).isFailure())
    return calDiagData; ///< return NULL ptr

  CalUtil::CalDiagnosticWord diagWord(trigBits, lacBits);

  /// reset return value to properly initialized CalDiagnosticData
  calDiagData.reset(new LdfEvent::CalDiagnosticData(diagWord.getDatum(), twr.val(), lyr.val()));

  return calDiagData;
}


/// calculate LAC bits for single layer
StatusCode CalDiagnosticTool::calcLACBits(const CalUtil::TwrNum twr,
                                          const CalUtil::LyrNum lyr,
                                          CalUtil::CalDiagnosticWord::CalDiagLACBits &lacBits) {

  /// compare signal in CIDAC vs THOLD for each crystal face (LRG
  /// diode only)
  for (FaceNum face; face.isValid(); face++)
    for (ColNum col; col.isValid(); col++) {
      const DiodeIdx diodeIdx(twr,lyr,col,face,LRG_DIODE);

      float signal;
      if (m_calSignalTool->getDiodeSignal(diodeIdx, signal).isFailure())
        return StatusCode::FAILURE;

      float thold;
      const FaceIdx faceIdx(twr,lyr,col,face);
      if (m_precalcCalibTool->getLacCIDAC(faceIdx, thold).isFailure())
        return StatusCode::FAILURE;

      if (signal > thold)
        lacBits[face][col] = true;
    }

  return StatusCode::SUCCESS;
}

/// calculate OR'd trigger bits for single cal diagnostic 
StatusCode CalDiagnosticTool::calcTrigBits(const CalUtil::TwrNum twr,
                                           const CalUtil::LyrNum lyr,
                                           CalUtil::CalDiagnosticWord::CalDiagTrigBits &trigBits) {

  for (FaceNum face; face.isValid(); face++) 
    for (DiodeNum diode; diode.isValid(); diode++)
      // 'or' all xtal faces on same GCRC together
      for (ColNum col; col.isValid(); col++) {
        bool trigBit;
        if (m_calTrigTool->getTriggerBit(DiodeIdx(twr,lyr,col,face,diode), trigBit).isFailure())
          return StatusCode::FAILURE;

        if (trigBit) {
          trigBits[face][diode] = true;
          break; // no need to process the rest of the xtals in this row
        }
      }

  return StatusCode::SUCCESS;
}

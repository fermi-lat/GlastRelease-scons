/**
 * @class CalClusterHitTool
 *
 * @brief Implements a Gaudi Tool for having access to the calorimeter hit
 *        xtals associated with a given cluster (based on the) information in
 *        the TDS relational tables.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"

#include "CalUtil/ICalClusterHitTool.h"


class CalClusterHitTool : public AlgTool, virtual public ICalClusterHitTool
{
public:
  /// Standard Gaudi Tool interface constructor
  CalClusterHitTool(const std::string& type,
                    const std::string& name,
                    const IInterface* parent);
  /// Destructor.
  virtual ~CalClusterHitTool() {}
  /// @brief Initialization of the tool.
  StatusCode initialize();
  /// @brief Fill the vectors of CalXtalRecData pointers for a given cluster.
  StatusCode fillRecDataVec(Event::CalCluster* cluster);
  /// @brief Return the m_recDataVec class member.
  inline std::vector<Event::CalXtalRecData*> getRecDataVec()
    const { return m_recDataVec; }
  /// @brief Set the m_recDataVec class member.
  void setRecDataVec(std::vector<Event::CalXtalRecData*> recDataVec)
    { m_recDataVec = recDataVec; }
  
private:

  /// Pointer to the data service.
  IDataProviderSvc* m_pEventSvc;
  /// Basic object of the tool (it stores the CalXtalRecData objects).
  std::vector<Event::CalXtalRecData*> m_recDataVec;
};


DECLARE_TOOL_FACTORY(CalClusterHitTool);


CalClusterHitTool::CalClusterHitTool(const std::string& type,
                                     const std::string& name,
                                     const IInterface* parent):
  AlgTool(type, name, parent)
{
  declareInterface<ICalClusterHitTool>(this);
  return;
}


StatusCode CalClusterHitTool::initialize()
{
  // Always believe in success (just kidding)!
  StatusCode sc = StatusCode::FAILURE;
  // Instantiate the message logger.
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "CalClusterHitTool is initializing..," << endreq;
  // Get a pointer to the event service.
  sc = serviceLocator()->service("EventDataSvc", m_pEventSvc, true);
  if (sc.isFailure()){
    log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
    return sc;
  }
  return sc;
}


StatusCode CalClusterHitTool::fillRecDataVec(Event::CalCluster* cluster)
{
  StatusCode sc = StatusCode::SUCCESS;
  // Clear the vector of CalXtalRecData pointers.
  m_recDataVec.clear();
  // Retrieve the list of relations (this is what we store in the TDS).
  Event::CalClusterHitTabList* tabList =
    SmartDataPtr<Event::CalClusterHitTabList>
    (m_pEventSvc, EventModel::CalRecon::CalClusterHitTab);
  // Construct the relation table from the list.
  Event::CalClusterHitTab* relTab = 0;
  if (tabList) {
    relTab = new Event::CalClusterHitTab(tabList);
  }
  // Finally, retrieve the list of associations and fill the m_recDataVec
  // data member.
  std::vector<Event::CalClusterHitRel*>relVec = relTab->getRelBySecond(cluster);
  if (!relVec.empty()) {
    std::vector<Event::CalClusterHitRel*>::const_iterator it;
    for (it = relVec.begin(); it != relVec.end(); it++){
      Event::CalXtalRecData* recData = (*it)->getFirst();
      m_recDataVec.push_back(recData);
    }
  }
  delete relTab;
  return StatusCode::SUCCESS;
}



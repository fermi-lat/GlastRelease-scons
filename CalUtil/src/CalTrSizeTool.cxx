/**
 * @class CalTrSizeTool
 *
 * @brief Implements a Gaudi Tool for evaluating the transverse size of
 *        a CAL cluster, give an input axis.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/TopLevel/EventModel.h"

#include "CalUtil/ICalTrSizeTool.h"


/// @brief Minimal utility class to store the basic xtal information
/// necessary to calculate the cumulative transverse size.

class CalTrSizeData
{
public:
  /// @brief Default constructor.
  CalTrSizeData(Event::CalXtalRecData* recData);    
  /// @brief Destructor.
  ~CalTrSizeData() {}
  /// @brief Set the reference axis for the object.
  ///
  /// If the transverse parameter is true only the transverse xtal information
  ///  is used.
  void setRefAxis(const Point& origin, const Vector& direction,
                  bool transverse = false);
  /// @brief Get the distance to the reference axis.
  inline double getDistToRefAxis() const { return m_distToRefAxis; }
  /// @brief Get the energy.
  inline double getEnergy() const { return m_energy; }
  /// @ Retireve the x, y, z coordinates of the hit.
  inline double x() const { return m_position.x(); }
  inline double y() const { return m_position.y(); }
  inline double z() const { return m_position.z(); }
  /// @brief Check wheter the underlying xtals lies along the x or yaxis.
  inline bool isx() const { return (m_layer % 2 == 0); }
  inline bool isy() const { return (m_layer % 2 != 0); }
  /// @brief Define how to sort.
  const bool operator<(const CalTrSizeData& other)  const
  { return m_distToRefAxis < other.getDistToRefAxis(); }
  /// @brief Standard output facility.
  std::ostream& fillStream(std::ostream& s) const;
  friend std::ostream& operator<< (std::ostream& s, const CalTrSizeData& obj)
  { return obj.fillStream(s); }
  
private:
  /// @brief hit position.
  Point m_position;
  /// @brief hit energy.
  double m_energy;
  /// @brief xtal tower.
  int m_tower;
  /// @brief xtal layer.
  int m_layer;
  /// @brief xtal column.
  int m_column;
  /// @brief The position of the xtal center.
  Point m_xtalCenter;
  /// @brief The direction of the xtal axis. 
  Vector m_xtalAxis;
  /// @brief Distance to the reference axis. 
  double m_distToRefAxis;
};


CalTrSizeData::CalTrSizeData(Event::CalXtalRecData* recData):
  m_position(recData->getPosition()),
  m_tower(recData->getTower()),
  m_layer(recData->getLayer()),
  m_column(recData->getColumn()),
  m_energy(recData->getEnergy()),
  m_distToRefAxis(0.)
{
  double towerPitch = 374.5;
  if (isx()) {
    m_xtalCenter = Point((m_tower - 4*(m_tower/4) - 1.5)*towerPitch, y(), z());
    m_xtalAxis = Vector(1., 0., 0.);
  } else {
    m_xtalCenter = Point(x(), (m_tower/4 - 1.5)*towerPitch, z());
    m_xtalAxis = Vector(0., 1., 0.);
  }
}


void CalTrSizeData::setRefAxis(const Point& origin, const Vector& direction,
                               bool transverse)
{
  if (transverse) {
    double halfXtalLength = 0.5*326.;
    // This is the slightly complicated case (credits Philippe Bruel). 
    Vector vec = m_xtalCenter - origin;
    double lambda = 1 - pow((m_xtalAxis*direction), 2.);
    if (lambda != 0) {
      lambda = (-(vec*m_xtalAxis) +
                (vec*direction)*(m_xtalAxis*direction))/lambda;
      if (lambda > halfXtalLength) lambda = halfXtalLength;
      else if (lambda < -halfXtalLength) lambda = -halfXtalLength;
      vec += lambda*m_xtalAxis;
    }
    double dist = vec.mag2() - pow((vec*direction), 2.);
    if (dist < 0.) dist = 0.;
    else dist = sqrt(dist);
    m_distToRefAxis = dist;
  } else {
    // And this is the easy case.
    m_distToRefAxis = (direction.cross(origin - m_position)).mag();
  }
}


std::ostream& CalTrSizeData::fillStream(std::ostream& s) const
{
  s << "Position = " << m_position << ", energy = " << m_energy <<
    ", d = " << m_distToRefAxis;
  return s;
}



/// @brief Another utility class to store the cumulative transverse size
/// distribution.

class CalCumulativeTrSize
{
public:
  /// @brief Default constructor.
  CalCumulativeTrSize() { clear(); }
  /// @brief Destructor.
  ~CalCumulativeTrSize() {}
  /// @brief Clear everything up.
  inline void clear()
  { m_fracVec.clear(); m_trSizeVec.clear(); }
  inline void addPoint(double frac, double trSize)
  { m_fracVec.push_back(frac); m_trSizeVec.push_back(trSize); }
  double getQuantile(double frac);
  
private:
  std::vector<double> m_fracVec;
  std::vector<double> m_trSizeVec;
};


double CalCumulativeTrSize::getQuantile(double frac)
{
  // Perform a linear interpolation of the underying data and return
  // the best estimate for the quantile.
  // Check if we're outside the range...
  if (frac <= m_fracVec.front()) {
    return m_trSizeVec.front();
  } else if (frac >= m_fracVec.back()) {
    return m_trSizeVec.back();
  } else {
    // ...and do a full binary search instead. 
    std::vector<double>::const_iterator lower =
      std::lower_bound(m_fracVec.begin(), m_fracVec.end(), frac);
    int index = int(lower - m_fracVec.begin());
    double x1 = m_fracVec.at(index - 1);
    double y1 = m_trSizeVec.at(index - 1);
    double x2 = m_fracVec.at(index);
    double y2 = m_trSizeVec.at(index);
    return y1 + (frac - x1)*(y2 - y1)/(x2 - x1);
  }
}



class CalTrSizeTool : public AlgTool, virtual public ICalTrSizeTool
{
public:
  /// Standard Gaudi Tool interface constructor
  CalTrSizeTool(const std::string& type,
                const std::string& name,
                const IInterface* parent);
  /// Destructor.
  virtual ~CalTrSizeTool() {}
  /// @brief Initialization of the tool.
  StatusCode initialize();
  /// @brief Fill the xtal collection.
  StatusCode fill(std::vector<Event::CalXtalRecData*> xtalList);
  StatusCode fill(Event::CalCluster* cluster);
  /// @brief Do the computation of the cumulative transverse size
  /// using the full xtal information.
  StatusCode computeTrSize(const Point& origin, const Vector& direction)
  { return compute(origin, direction, false); }
  /// @brief Do the computation of the cumulative transverse size
  /// using only the transverse xtal information.
  StatusCode computeTrSizeTrans(const Point& origin, const Vector& direction)
  { return compute(origin, direction, true); }
  /// @brief Get a generic quantile of the cumulative transverse size.
  double getQuantile(double frac)
  { return m_cumulativeTrSize.getQuantile(frac); }

private:
  /// @brief Pointer to the data service.
  IDataProviderSvc* m_pEventSvc;
  /// @brief The underlying CalTrSizeData objects.
  std::vector<CalTrSizeData> m_trSizeData;
  /// @brief The object storing the cumulative transverse size.
  CalCumulativeTrSize m_cumulativeTrSize;
  /// @brief The total energy in the xtal collection.
  double m_totalEnergy;
  /// @brief Do the actual computation of the cumulative transverse size.
  StatusCode compute(const Point& origin, const Vector& direction,
                     bool transverse);
};


DECLARE_TOOL_FACTORY(CalTrSizeTool);


CalTrSizeTool::CalTrSizeTool(const std::string& type,
                             const std::string& name,
                             const IInterface* parent):
  AlgTool(type, name, parent)
{
  declareInterface<ICalTrSizeTool>(this);
  return;
}


StatusCode CalTrSizeTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  // Instantiate the message logger.
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "CalTrSizeTool is initializing..," << endreq;
  // Get a pointer to the event service.
  sc = serviceLocator()->service("EventDataSvc", m_pEventSvc, true);
  if (sc.isFailure()){
    log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
  }
  return sc;
}


StatusCode CalTrSizeTool::fill(std::vector<Event::CalXtalRecData*> xtalList)
{
  // Clear the underlying xtal vector...
  m_trSizeData.clear();
  m_totalEnergy = 0.;
  if (xtalList.size() == 0) return StatusCode::FAILURE;
  // ...and fill it starting from the xtal collection.
  std::vector<Event::CalXtalRecData*>::const_iterator xtal;
  for (xtal = xtalList.begin(); xtal != xtalList.end(); xtal++) {
    m_trSizeData.push_back(CalTrSizeData(*xtal));
    m_totalEnergy += (*xtal)->getEnergy();
  }
  return StatusCode::SUCCESS;
}


StatusCode CalTrSizeTool::fill(Event::CalCluster* cluster)
{
  StatusCode sc = StatusCode::SUCCESS;
  std::vector<Event::CalXtalRecData*> xtalList;
  xtalList.clear();
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
      xtalList.push_back(recData);
    }
  }
  delete relTab;
  return fill(xtalList);
}


StatusCode CalTrSizeTool::compute(const Point& origin, const Vector& direction,
                                  bool transverse)
{
  m_cumulativeTrSize.clear();
  double runningEnergy = 0.;
  double runningTrSize = 0.;
  if (m_trSizeData.size() == 0) return StatusCode::FAILURE;
  // Set the reference axis for the underlying CalTrSize data obejcts.
  for (std::vector<CalTrSizeData>::iterator xtal = m_trSizeData.begin();
       xtal != m_trSizeData.end(); xtal++) {
    (*xtal).setRefAxis(origin, direction, transverse);
  }
  // Sort the xtals according to their distance from the reference axis.
  std::sort(m_trSizeData.begin(), m_trSizeData.end());
  // Loop again to build the integral curve.
  for (std::vector<CalTrSizeData>::const_iterator xtal = m_trSizeData.begin();
       xtal != m_trSizeData.end(); xtal++) {
    //std::cout << *xtal << std::endl;
    double energy = (*xtal).getEnergy();
    double dist = (*xtal).getDistToRefAxis();
    runningEnergy += energy;
    runningTrSize += dist*dist*energy;
    double frac = runningEnergy/m_totalEnergy;
    double trSize = sqrt(runningTrSize/runningEnergy);
    m_cumulativeTrSize.addPoint(frac, trSize);
  }
  return StatusCode::SUCCESS;
}

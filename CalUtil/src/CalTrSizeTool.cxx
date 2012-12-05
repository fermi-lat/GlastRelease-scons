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
  void setRefAxis(const Point& origin, const Vector& direction);
  /// @brief Get the distance to the reference axis.
  inline double getDistToRefAxis() const { return m_distToRefAxis; }
  /// @brief Get the energy.
  inline double getEnergy() const { return m_energy; }
  /// @brief Define how to sort.
  const bool operator<(const CalTrSizeData& other)  const
  { return m_distToRefAxis < other.getDistToRefAxis(); }
  /// @brief Standard output facility.
  std::ostream& fillStream(std::ostream& s) const;
  friend std::ostream& operator<< (std::ostream& s, const CalTrSizeData& obj)
  { return obj.fillStream(s); }
  
private:
  /// @brief xtal position.
  Point m_position;
  /// @brief xtal energy.
  double m_energy;
  /// @brief Distance to the reference axis. 
  double m_distToRefAxis;
};


CalTrSizeData::CalTrSizeData(Event::CalXtalRecData* recData):
  m_position(recData->getPosition()),
  m_energy(recData->getEnergy()),
  m_distToRefAxis(0.)
{}


void CalTrSizeData::setRefAxis(const Point& origin, const Vector& direction)
{
  m_distToRefAxis = (direction.cross(origin - m_position)).mag();
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
    //return 0.5*(m_trSizeVec.at(index-1) + m_trSizeVec.at(index));
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
  /// @brief 
  StatusCode fill(std::vector<Event::CalXtalRecData*> xtalList);
  /// @brief
  StatusCode calculate(const Point& origin, const Vector& direction);
  /// @brief 
  double getQuantile(double frac)
  { return m_cumulativeTrSize.getQuantile(frac); }

private:
  /// @brief The underlying CalTrSizeData objects.
  std::vector<CalTrSizeData> m_trSizeData;
  /// @brief The object storing the cumulative transverse size.
  CalCumulativeTrSize m_cumulativeTrSize;
  /// @brief The total energy in the xtal collection.
  double m_totalEnergy;
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
  return sc;
}


StatusCode CalTrSizeTool::fill(std::vector<Event::CalXtalRecData*> xtalList)
{
  // Clear the underlying xtal vector...
  m_trSizeData.clear();
  m_totalEnergy = 0.;
  // ...and fill it starting from the xtal collection.
  std::vector<Event::CalXtalRecData*>::const_iterator xtal;
  for (xtal = xtalList.begin(); xtal != xtalList.end(); xtal++) {
    m_trSizeData.push_back(CalTrSizeData(*xtal));
    m_totalEnergy += (*xtal)->getEnergy();
  }
  return StatusCode::SUCCESS;
}


StatusCode CalTrSizeTool::calculate(const Point& origin,
                                    const Vector& direction)
{
  m_cumulativeTrSize.clear();
  double runningEnergy = 0.;
  double runningTrSize = 0.;
  // Set the reference axis for the underlying CalTrSize data obejcts.
  for (std::vector<CalTrSizeData>::iterator xtal = m_trSizeData.begin();
       xtal != m_trSizeData.end(); xtal++) {
    (*xtal).setRefAxis(origin, direction);
  }
  // Sort the xtals according to their distance from the reference axis.
  std::sort(m_trSizeData.begin(), m_trSizeData.end());
  // Loop again to build the integral curve.
  for (std::vector<CalTrSizeData>::const_iterator xtal = m_trSizeData.begin();
       xtal != m_trSizeData.end(); xtal++) {
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

/**
 * @class CalDirFittingTool
 *
 * @brief Implements a Gaudi Tool for fitting the CAL direction in
 *        several different flavours.
 *
 * @author Luca Baldini (luca.baldini@pi.infn.it)
 *
 * $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/GaudiException.h"

#include "CalUtil/ICalDirFittingTool.h"


/// Convenience function to check whether a given layer is x or y.
/// A "x" xtal *lies* along the x direction, i.e. the transverse coordinate is y.
inline bool isx(int layerId) { return layerId % 2 == 0; }
inline bool isy(int layerId) { return layerId % 2 != 0; }


class CalDirFittingTool : public AlgTool, virtual public ICalDirFittingTool
{
public:
  /// Standard Gaudi Tool interface constructor
  CalDirFittingTool(const std::string& type,
                    const std::string& name,
                    const IInterface* parent);
  /// Destructor.
  virtual ~CalDirFittingTool() {}
  /// @brief Initialization of the tool.
  StatusCode initialize();
  /// @brief Fit using only the transverse information separately on the
  ///        x-z and y-z projections.
  StatusCode transverseFit2d(Event::CalCluster* cluster, double powerWeight = 1.);
  inline Event::CalFitParams getFitParams() const { return m_fitParams; }

private:
  /// @brief The underlying CalFitParams object to store the fit results.
  Event::CalFitParams m_fitParams;
  /// @brief Convert the fit parameters from the slope/intercept representation
  ///        to the cosine directions and fill the m_fitParams class member.
  void setFitParams(double xzintercept, double xzslope,
                    double yzintercept, double yzslope,
                    double zcentroid);
};


DECLARE_TOOL_FACTORY(CalDirFittingTool);


CalDirFittingTool::CalDirFittingTool(const std::string& type,
                                     const std::string& name,
                                     const IInterface* parent):
  AlgTool(type, name, parent)
{
  declareInterface<ICalDirFittingTool>(this);
  return;
}


StatusCode CalDirFittingTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;
  // Instantiate the message logger.
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "CalDirFittingTool is initializing..," << endreq;
  return sc;
}


void CalDirFittingTool::setFitParams(double xzintercept, double xzslope,
                                     double yzintercept, double yzslope,
                                     double zcentroid)
{
  /// To be implemented.
}


StatusCode CalDirFittingTool::transverseFit2d(Event::CalCluster* cluster,
                                              double powerWeight)
{
  // Clear the CalFitParams container.
  m_fitParams.clear();
  
  // Initialize some variables.
  int nx, ny = 0;
  double x0, y0, z0, wx, wy, wz = 0.;
  double sx, sx_z, sx_zz, sx_x, sx_zx = 0.;
  double sy, sy_z, sy_zz, sy_y, sy_zy = 0.;
  double totalEnergy = (*cluster).getXtalsParams().getXtalRawEneSum();

  // Loop over the layers.
  for (int layerId = 0; layerId < NUMCALLAYERS; layerId++) {
    double energy = (*cluster)[layerId].getEnergy();
    if (energy > 0) {
      Point position = (*cluster)[layerId].getPosition();
      double err2 =  pow((totalEnergy*energy), powerWeight);
      if (isx(layerId)){
        ny += 1;
        y0 += energy*position.y();
        wy += energy;
        sy += 1./err2;
        sy_z += position.z()/err2;
        sy_zz += position.z()*position.z()/err2;
        sy_y += position.y()/err2;
        sy_zy += position.z()*position.y()/err2;
      } else {
        nx += 1;
        x0 += energy*position.x();
        wx += energy;
        sx += 1./err2;
        sx_z += position.z()/err2;
        sx_zz += position.z()*position.z()/err2;
        sx_x += position.x()/err2;
        sx_zx += position.z()*position.x()/err2;
      }
      z0 += energy*position.z();
      wz += energy;
    }
  }

  // Need at least two points in both views.
  if (nx < 2 || ny < 2) return StatusCode::FAILURE;

  // Normalize the components of the centroid position.
  if (wx > 0) x0 /= wx;
  if (wy > 0) y0 /= wy;
  if (wz > 0) z0 /= wz;
  // Go ahead with the fit parameters: these are the slopes and intercepts
  // in the two views.
  double xzintercept = (sx_x*sx_zz - sx_zx*sx_z)/(sx_zz*sx - sx_z*sx_z);
  double xzslope = (sx_zx*sx - sx_z*sx_x)/(sx_zz*sx - sx_z*sx_z);
  double yzintercept = (sy_y*sy_zz - sy_zy*sy_z)/(sy_zz*sy - sy_z*sy_z);
  double yzslope = (sy_zy*sy - sy_z*sy_y)/(sy_zz*sy - sy_z*sy_z);
  // Now we want to refer the centroid to the common z0 value
  // calculated in the previous loop.
  double xcen = xzintercept + z0*xzslope;
  double ycen = yzintercept + z0*yzslope; 
  double zcen = z0;
  Point centroid(xcen, ycen, zcen);
  // Finally: convert from the slope/intercept representation to the
  // director cosines.
  double denom = sqrt(1 + xzslope*xzslope + yzslope*yzslope);
  double xdir  = xzslope/denom;
  double ydir  = yzslope/denom;
  double zdir  = 1./denom;
  Vector axis(xdir, ydir, zdir);

  // One more loop to calculate the chisquare.
  double chisq = 0.;
  for (int layerId = 0; layerId < NUMCALLAYERS; layerId++) {
    double energy = (*cluster)[layerId].getEnergy();
    if (energy > 0) {
      Point position = (*cluster)[layerId].getPosition();
      double err2 =  pow((totalEnergy*energy), powerWeight);
      if (isx(layerId)){
        double dy = position.y() - (ycen + (position.z() - zcen)/ydir);
        chisq += dy*dy/err2;
      } else {
        double dx = position.x() - (xcen + (position.z() - zcen)/xdir);;
        chisq += dx*dx/err2;
      }
    }
  }

  // Update the CalFitParams container.
  m_fitParams.setCentroid(centroid);
  m_fitParams.setAxis(axis);
  m_fitParams.setFitLayers(nx + ny);
  m_fitParams.setChiSquare(chisq);

  return StatusCode::SUCCESS;
}

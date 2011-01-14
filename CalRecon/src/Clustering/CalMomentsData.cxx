/**
   @file CalMomentsData.h
   
   @brief Implementation file for the CalMomentsData class.

   @author Luca Baldini (luca.baldini@pi.infn.it).
   
   Revision $Revision$, commited on $Date$.
   $Id$
*/

#include "src/Clustering/CalMomentsData.h"

#include <algorithm>


CalMomentsData::CalMomentsData(const Point& position, const double weight,
			       int tower, int layer, int column) :
  m_statusBits(ZERO),
  m_basePosition(position),
  m_corrPosition(Point(-9999., -9999., -9999.)),
  m_weight(weight),
  m_tower(tower),
  m_layer(layer),
  m_column(column),
  m_distToAxis(0.),
  m_coordAlongAxis(0.)
{ 
  // Nothing to do, here.
}

CalMomentsData::CalMomentsData(Event::CalXtalRecData* recData) :
  m_statusBits(ZERO),
  m_corrPosition(Point(-9999., -9999., -9999.)),
  m_distToAxis(0.),
  m_coordAlongAxis(0.)
{
  m_basePosition = recData->getPosition();
  m_weight = recData->getEnergy();
  m_tower = (recData->getPackedId()).getTower();
  m_layer = (recData->getPackedId()).getLayer();
  m_column = (recData->getPackedId()).getColumn();
}

const Point& CalMomentsData::getPosition() const
{
  if ( checkStatusBit(USE_FIT_POS) ) {
    return m_corrPosition;
  }
  else {
    return m_basePosition;
  }
}

double CalMomentsData::getFitCorrAmount() const
{
  if ( checkStatusBit(FIT_POS_AVAILABLE) ) {
    if ( isx() ) {
      return m_corrPosition.x() - m_basePosition.x();
    }
    else {
      return m_corrPosition.y() - m_basePosition.y();
    }
  }
  else {
    return -9999.;
  }
}

void CalMomentsData::applyFitCorrection(const Event::CalFitParams fitParams,
					const ICalReconSvc* calReconSvc)
{
  Point fitCentroid = fitParams.getCentroid();
  Vector fitAxis = fitParams.getAxis();
  int numFitLayers = fitParams.getFitLayers();

  // If there's no useful fit information, give up on the correction.
  if ( (numFitLayers == 0 ) || (fitAxis.z() == 0.) ) {
    return;
  }
  
  // Otherwise carry on bravely and store the relevant information.
  double towerPitch = calReconSvc->getCaltowerPitch();
  double csILength  = calReconSvc->getCalCsILength();
  bool flightGeom   = calReconSvc->getCalFlightGeom();
  double x = m_basePosition.x();
  double y = m_basePosition.y();
  double z = m_basePosition.z();
  int towerIdy = m_tower / 4;
  int towerIdx = m_tower - 4*towerIdy;
  double minCoord;
  double maxCoord;
  double distToEdge;
  // Xtals along the x coordinate first...
  if ( isx() ) {
    minCoord = -1.5*towerPitch + towerPitch*(double)towerIdx - csILength/2;
    maxCoord = -1.5*towerPitch + towerPitch*(double)towerIdx + csILength/2;
    distToEdge = std::min( fabs(x - minCoord), fabs(x - maxCoord) );
    if ( distToEdge < DIST_LONG_POS_INVALID ) {
      setStatusBit(LONG_POS_INVALID);
    }
    x = fitCentroid.x() + (z - fitCentroid.z()) / fitAxis.z() * fitAxis.x();
    distToEdge = std::min( fabs(x - minCoord), fabs(x - maxCoord) );
    if ( (x > minCoord) && (x < maxCoord) && (distToEdge < DIST_FIT_POS_NEAR_EDGE) ) {
      setStatusBit(FIT_POS_NEAR_EDGE);
    }
    if ( x < minCoord ) {
      x = minCoord;
      setStatusBit(FIT_POS_INVALID);
    }
    if ( x > maxCoord ) {
      x = maxCoord;
      setStatusBit(FIT_POS_INVALID);
    }
    m_corrPosition = Point(x, y, z);
  }
  // ... and then xtals along the y coordinate.
  else {
    // And this is slightly more complicated as we have to distinguish between the LAT...
    if ( flightGeom ) {
      minCoord = -1.5*towerPitch + towerPitch*(double)towerIdy - csILength/2;
      maxCoord = -1.5*towerPitch + towerPitch*(double)towerIdy + csILength/2;
    }
    // ...and the CU.
    else {
      minCoord = - csILength/2;
      maxCoord = csILength/2;
    }
    distToEdge = std::min( fabs(y - minCoord), fabs(y - maxCoord) );
    if ( distToEdge < DIST_LONG_POS_INVALID ) {
      setStatusBit(LONG_POS_INVALID);
    }
    y = fitCentroid.y() + (z - fitCentroid.z()) / fitAxis.z() * fitAxis.y();
    distToEdge = std::min( fabs(y - minCoord), fabs(y - maxCoord) );
    if ( (y > minCoord) && (y < maxCoord) && (distToEdge < DIST_FIT_POS_NEAR_EDGE) ) {
      setStatusBit(FIT_POS_NEAR_EDGE);
    }
    if ( y < minCoord ) {
      y = minCoord;
      setStatusBit(FIT_POS_INVALID);
    }
    if ( y > maxCoord ) {
      y = maxCoord;
      setStatusBit(FIT_POS_INVALID);
    }
    m_corrPosition = Point(x, y, z);
  }
  setStatusBit(FIT_POS_AVAILABLE);
}

void CalMomentsData::enableFitCorrection()  
{
  if ( checkStatusBit(FIT_POS_AVAILABLE) && !checkStatusBit(FIT_POS_INVALID) ) {
    setStatusBit(USE_FIT_POS);
  }
}

double CalMomentsData::calcDistToAxis(const Point& centroid, const Vector& axis)
{
  Vector diffVec   = centroid - getPosition();
  Vector crossProd = axis.cross(diffVec);
  return m_distToAxis = crossProd.mag();
}

double CalMomentsData::calcCoordAlongAxis(const Point& centroid, const Vector& axis)
{
  return m_coordAlongAxis = (getPosition() - centroid).dot(axis);
}

void CalMomentsData::setReferenceAxis(const Point& centroid, const Vector& axis)
{
  m_distToAxis = axis.cross(centroid - getPosition()).mag();
  m_coordAlongAxis = (getPosition() - centroid).dot(axis);
}

std::ostream& CalMomentsData::fillStream(std::ostream& s) const
{
  s <<
    "Status bits: 0x" << std::hex << m_statusBits << std::dec << "\n"
    "Tower = " << m_tower << ", layer = " << m_layer << ", column = " << m_column << "\n" <<
    "Position = " << m_basePosition << ", corrected position = " << m_corrPosition <<
    ", correction = " << getFitCorrAmount();
  return s;
}

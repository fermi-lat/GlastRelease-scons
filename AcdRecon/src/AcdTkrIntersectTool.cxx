#include "AcdTkrIntersectTool.h"
#include "GaudiKernel/MsgStream.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

#include "Event/Recon/TkrRecon/TkrTrackHit.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

#include "idents/AcdId.h"

DECLARE_TOOL_FACTORY(AcdTkrIntersectTool)

AcdTkrIntersectTool::AcdTkrIntersectTool
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent )
 { declareInterface<AcdITkrIntersectTool>(this) ; }

AcdTkrIntersectTool::~AcdTkrIntersectTool()
{} 

// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode AcdTkrIntersectTool::initialize()
 {
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;

  // GlastDetSvc
  sc = service("GlastDetSvc", m_detSvc);
  if (sc.isFailure()) {
    log<< MSG::ERROR<<"GlastDetSvc not found"<<endreq ;
    return StatusCode::FAILURE ;
  } 
  
  // TkrGeometrySvc
  sc = service("TkrGeometrySvc", m_tkrGeom, true);
  if ( sc.isFailure()) {
    log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
    return StatusCode::FAILURE;
  }

  // G4Propagator
  if( ! toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
    log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
    return StatusCode::FAILURE;
  }
  
  log<<MSG::INFO<<"END initialize()"<<endreq ;
  return StatusCode::SUCCESS ;
 }

StatusCode AcdTkrIntersectTool::findIntersections (const Event::TkrTrackCol * tracks,
						   Event::AcdTkrIntersectionCol * intersects,
						   std::map<idents::AcdId,unsigned char>& hitMap)
  //Purpose and method:
  //
  // TDS input: TkrTrackCol
  // TDS output: AcdTkrIntersectionCol
{
  MsgStream log(msgSvc(),name()) ;
  output = intersects;

  // loop on tracks
  int iTrack(0);
  for ( Event::TkrTrackCol::const_iterator it = tracks->begin();
	it != tracks->end(); it++, iTrack++ ) {
    const Event::TkrTrack* aTrack = *it;
    if ( aTrack == 0 ) return StatusCode::FAILURE ;
    // do the intersections for this track
    int numberIntersections = doTrack(*aTrack,iTrack,hitMap,log);    
    if ( numberIntersections < 0 ) {
      // negative return values indicate failure
      return StatusCode::FAILURE;
    }
  }
  // done
  return StatusCode::SUCCESS ;
}

int AcdTkrIntersectTool::doTrack(const Event::TkrTrack& aTrack, int iTrack, 
				 std::map<idents::AcdId,unsigned char>& hitMap, MsgStream& log) {
  
  // Define the fiducial volume of the LAT
  // FIXME -- this should come for some xml reading service
  //
  // top is defined by planes at + 754.6 -> up to stacking of tiles
  // sides are defined by planes at +-840.14
  // add 10 cm to make sure that we catch everything
  static double top_distance = 854.6;     // center of tiles in cols 1 and 3
  static double side_distance = 940.14;   // center of tiles in sides

  // we will return the number of intersection w/ the Acd for this track
  int returnValue(0);

  // Get the start point & direction of the track & the params & energy also
  const Point initialPosition = aTrack.getInitialPosition();
  const Vector initialDirection = aTrack.getInitialDirection();
  const Event::TkrTrackParams& trackPars = aTrack[0]->getTrackParams(Event::TkrTrackHit::SMOOTHED);
  double startEnergy = aTrack[0]->getEnergy();

  // hits -x or +x side ?
  const double normToXIntersection =  initialDirection.x() > 0 ?  
    -1.*side_distance - initialPosition.x() :    // hits -x side
    1.*side_distance - initialPosition.x();      // hits +x side  
  const double slopeToXIntersection = fabs(initialDirection.x()) > 1e-9 ? 
    -1. / initialDirection.x() : 1e9;
  const double sToXIntersection = normToXIntersection * slopeToXIntersection;

  // hits -y or +y side ?
  const double normToYIntersection = initialDirection.y() > 0 ?  
    -1.*side_distance - initialPosition.y() :    // hits -y side
    1.*side_distance - initialPosition.y();      // hits +y side
  const double slopeToYIntersection = fabs(initialDirection.y()) > 1e-9 ? 
    -1. / initialDirection.y() : 1e9;
  const double sToYIntersection = normToYIntersection * slopeToYIntersection;

  // hits top
  const double normToZIntersection = top_distance - initialPosition.z();
  const double slopeToZIntersection = -1. / initialDirection.z();
  const double sToZIntersection = normToZIntersection * slopeToZIntersection;

  // pick closest plane
  const double arcLength = sToXIntersection < sToYIntersection ?
    ( sToXIntersection < sToZIntersection ? sToXIntersection : sToZIntersection ) :
    ( sToYIntersection < sToZIntersection ? sToYIntersection : sToZIntersection ) ;
  
  // protect against negative arcLenghts
  if ( arcLength < 0. ) return -1;

  // setup the G4Propagator
  bool downwards = initialDirection.z() < 0.0;
  m_G4PropTool->setStepStart(trackPars,initialPosition.z(),downwards); 
  m_G4PropTool->step(arcLength);  

  // setup some of the variables for the step loop
  idents::VolumeIdentifier volId;
  idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
  double arcLengthToIntersection(0.);

  // do the step loop
  int numSteps = m_G4PropTool->getNumberSteps();  
  for(int istep = 0; istep < numSteps; ++istep) {
    // running sums of arcLength
    double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
    arcLengthToIntersection += arcLen_step;
    // check to see if this step is in the Acd
    volId = m_G4PropTool->getStepVolumeId(istep);
    volId.prepend(prefix);
    if ( ! volId.isAcd() ) continue;
    switch ( volId[2] ) {
    case 40: // Acd Tile
    case 41: // Acd ribbon
      if ( volId.size() < 5 ) {
	// FIXME -- need to verify that this is really dead material
	continue;
      }
      break;
    default:
      // Dead Acd material
      continue;
    }
    idents::AcdId acdId(volId);

    // check if this tile was hit or not
    unsigned char hitMask(0);    
    std::map<idents::AcdId,unsigned char>::const_iterator findHit = hitMap.find(acdId);
    if ( findHit != hitMap.end() ) { hitMask = findHit->second; }

    // which face of detector (0=top,1=-x,2=-y,3=+x,4=+y)
    int face = volId[1];

    // query the detector service for the local frame
    HepTransform3D transform; 
    StatusCode sc = m_detSvc->getTransform3DByID(volId, &transform);
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Could not get transform for id " << volId.name() << endreq;
      return -1;
    }
    HepPoint3D center(0., 0., 0.);
    HepPoint3D xT = transform * center;
    
    double localPosition[2];
    HepMatrix localCovMatrix(2,2);
    double check(0.);  // this should usually be 1/2 the tile thickness
    double delta(0.);

    // propgate the postion and error matrix to the intersection
    Point x_step       = m_G4PropTool->getStepPosition(istep);
    Event::TkrTrackParams next_params = m_G4PropTool->getTrackParams(arcLengthToIntersection,startEnergy,true);

    // which face of the ACD was hit?
    // calcluate the local variable correctly for that face
    switch ( face ) {
    case 0:
      // top: x,y,z are same for local and global
      localPosition[0] = x_step.x() - xT.x();
      localPosition[1] = x_step.y() - xT.y();
      check = x_step.z() - xT.z();
      delta = xT.z() - initialPosition.z();
      errorAtZPlane(delta,next_params,localCovMatrix);
      break;
    case 1:
    case 3:
      // +-x sides: y -> local x, z -> local y
      localPosition[0] = x_step.y() - xT.y();
      localPosition[1] = x_step.z() - xT.z();
      check = x_step.x() - xT.x();
      delta = xT.x() - initialPosition.x();
      errorAtXPlane(delta,next_params,localCovMatrix);
      break;
    case 2:
    case 4:
      // +-y sides: x -> local x, z -> local y
      localPosition[0] = x_step.x() - xT.x();
      localPosition[1] = x_step.z() - xT.z();
      check = x_step.y() - xT.y();	    
      delta = xT.y() - initialPosition.y();
      errorAtYPlane(delta,next_params,localCovMatrix);
      break;
    } 

    if ( false ) {
      double localXErr = sqrt(localCovMatrix[0][0]);
      double localYErr = sqrt(localCovMatrix[1][1]);
      double correl =  localCovMatrix[0][1] / ( localXErr * localYErr );    
      
      std::cout << istep << ' ' << volId.name() << ' ' << acdId.id() << ' ' << acdId.volId().name() << ' '
		<< arcLen_step 
		<< " (" << x_step.x() << ',' << x_step.y() << ',' << x_step.z() 
		<< ") [" << localPosition[0] << ',' << localPosition[1] << ',' << check << "] " 
		<< delta << ' ' << arcLengthToIntersection << ' ' << startEnergy
		<< " {" << localXErr << ',' << localYErr << ',' << correl << "}. "
		<< hitMask
		<< std::endl;
    }      

    // ok, we have everything we need, build the AcdTkrIntersection 
    Event::AcdTkrIntersection* theIntersection = 
      new Event::AcdTkrIntersection(acdId,iTrack,
				    x_step,
				    localPosition,localCovMatrix,
				    arcLengthToIntersection,arcLen_step,
				    hitMask);
    // print to the log if debug level is set
    theIntersection->writeOut(log);
    
    // append to collection and increment counter
    output->push_back(theIntersection);
    returnValue++;
  }

  return returnValue;
}

void AcdTkrIntersectTool::errorAtXPlane(const double delta, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) const {
  
  // get the tk params.  
  // want the x and y and the normal to the plane
  const double m_x = pars.getxSlope();
  const double inv_m_x = fabs(m_x) > 1e-9 ? 1./ m_x : 1e9;
  const double m_y = pars.getySlope();

  // ok input the jacobian
  //
  //     The intersection occurs at:
  //       I_y = y_0 - delta * m_y / m_x
  //       I_z = z_0 - delta * 1. / m_x
  //     where: 
  //       delta = xPlane - x_0
  //
  //     The jacobian terms are defined as:
  //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
  //  
  HepMatrix jacobian(2,4);
  //jacobian[0][0] = 0.;                           // dI_y / dx_0  = 0
  //jacobian[0][1] = -delta * inv_m_x;             // dI_y / dm_x  = -d / m_x
  //jacobian[0][2] = 1.;                           // dI_y / dy_0  = 1
  //jacobian[0][3] = delta * m_y * inv_m_x;        // dI_y / dm_y  = -d * m_y / m_x

  //jacobian[1][0] = - inv_m_x;                    // dI_z / dx_0  = -1 / m_x
  //jacobian[1][1] = delta * inv_m_x * inv_m_x;    // dI_z / dm_x  = -d / m_x * m_x
  //jacobian[1][2] = 0.;                           // dI_z / dy_0  = 0
  //jacobian[1][3] = 0.;                           // dI_z / dm_y  = 0

  jacobian[0][0] = 0.;                           // dI_y / dx_0  = 0
  jacobian[0][1] = 0.;
  jacobian[0][2] = 1.;                           // dI_y / dy_0  = 1
  jacobian[0][3] = 0.;

  jacobian[1][0] = - inv_m_x;                    // dI_z / dx_0  = -1 / m_x
  jacobian[1][1] = delta * inv_m_x * inv_m_x;    // dI_z / dm_x  = -d / m_x * m_x
  jacobian[1][2] = 0.;                           // dI_z / dy_0  = 0
  jacobian[1][3] = 0.;                           // dI_z / dm_y  = 0  

  for ( unsigned int i(0);  i < 2; i++ ) {
    for ( unsigned int j(0);  j < 2; j++ ) {
      double currentSum = 0.;
      for ( unsigned int k(1);  k < 5; k++ ) {
	for ( unsigned int l(1);  l < 5; l++ ) {
	  currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	}
      }
      covAtPlane[i][j] = currentSum;
    }
  }
}

void AcdTkrIntersectTool::errorAtYPlane(const double delta, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) const {

  // get the tk params.  
  // want the x and y and the normal to the plane
  const double m_x = pars.getxSlope();
  const double m_y = pars.getySlope();
  const double inv_m_y = fabs(m_y) > 1e-9 ? 1./ m_y : 1e9;

  // ok input the jacobian
  //
  //     The intersection occurs at:
  //       I_x = x_0 - delta * m_x / m_y
  //       I_z = z_0 - delta * 1. / m_y
  //     where:
  //       delta = yPlane - y_0
  //
  //     The jacobian terms are defined as:
  //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
  //  
  HepMatrix jacobian(2,4);
  //jacobian[0][0] = 1.;                           // dI_x / dx_0  = 1                        
  //jacobian[0][1] = delta * m_x * inv_m_y;        // dI_x / dm_x  = d m_x / m_y
  //jacobian[0][2] = 0.;                           // dI_x / dy_0  = 0
  //jacobian[0][3] = -delta * inv_m_y;             // dI_x / dm_y  = -d / m_y

  //jacobian[1][0] = 0.;                           // dI_z / dx_0  = 0                    
  //jacobian[1][1] = 0.;                           // dI_z / dm_x  = 0                    
  //jacobian[1][2] = - inv_m_y;                    // dI_z / dy_0  = -1 / m_y
  //jacobian[1][3] = delta * inv_m_y * inv_m_y;    // dI_z / dm_y  = d / m_y * m_y
  
  jacobian[0][0] = 1.;                           // dI_x / dx_0  = 1                        
  jacobian[0][1] = 0.;
  jacobian[0][2] = 0.;                           // dI_x / dy_0  = 0
  jacobian[0][3] = 0.;

  jacobian[1][0] = 0.;                           // dI_z / dx_0  = 0                    
  jacobian[1][1] = 0.;                           // dI_z / dm_x  = 0                    
  jacobian[1][2] = - inv_m_y;                    // dI_z / dy_0  = -1 / m_y
  jacobian[1][3] = delta * inv_m_y * inv_m_y;    // dI_z / dm_y  = d / m_y * m_y  

  for ( unsigned int i(0);  i < 2; i++ ) {
    for ( unsigned int j(0);  j < 2; j++ ) {
      double currentSum = 0.;
      for ( unsigned int k(1);  k < 5; k++ ) {
	for ( unsigned int l(1);  l < 5; l++ ) {
	  currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	}
      }
      covAtPlane[i][j] = currentSum;
    }
  }
}

void AcdTkrIntersectTool::errorAtZPlane(const double delta, const Event::TkrTrackParams& pars, HepMatrix& covAtPlane) const {

  // ok input the jacobian
  //
  //     The intersection occurs at:
  //       I_x = x_0 - delta * m_x
  //       I_y = y_0 - delta * m_y
  //     where:
  //       delta = zPlane - z_0
  //
  //     The jacobian terms are defined as:
  //        J[i][k] = d I_i / d param_k     where param = { x, y, m_x, m_y }
  //  
  HepMatrix jacobian(2,4);
  
  //jacobian[0][0] = 1.;            // dI_x / dx_0 = 1
  //jacobian[0][1] = -delta;        // dI_x / dm_x = -d
  //jacobian[0][2] = 0.;            // dI_x / dy_0 = 0
  //jacobian[0][3] = 0.;            // dI_x / dm_y = 0

  //jacobian[1][0] = 0.;            // dI_y / dx_0 = 0
  //jacobian[1][1] = 0.;            // dI_y / dm_x = 0
  //jacobian[1][2] = 1.;            // dI_y / dy_0 = 1
  //jacobian[1][3] = -delta;        // dI_y / dm_y = -d

  jacobian[0][0] = 1.;            // dI_x / dx_0 = 1
  jacobian[0][1] = 0.;        // dI_x / dm_x = -d
  jacobian[0][2] = 0.;            // dI_x / dy_0 = 0
  jacobian[0][3] = 0.;            // dI_x / dm_y = 0

  jacobian[1][0] = 0.;            // dI_y / dx_0 = 0
  jacobian[1][1] = 0.;            // dI_y / dm_x = 0
  jacobian[1][2] = 1.;            // dI_y / dy_0 = 1
  jacobian[1][3] = 0.;        // dI_y / dm_y = -d  

  for ( unsigned int i(0);  i < 2; i++ ) {
    for ( unsigned int j(0);  j < 2; j++ ) {
      double currentSum = 0.;
      for ( unsigned int k(1);  k < 5; k++ ) {
	for ( unsigned int l(1);  l < 5; l++ ) {
	  currentSum += jacobian[i][k-1] * jacobian[j][l-1] * pars(k,l);
	}
      }
      covAtPlane[i][j] = currentSum;
    }
  }
}

/** 
* @class G4PropagationTool
*
* @brief A Gaudi Tool for propagating a track through the GLAST detector. 
*
*        This tool implements a "particle propagator" for transporting a track
*        through the GLAST detector. It inherits from the class ParticleTransporter,
*        which does the actual tracking using Geant4 volume routines and the detmodel
*        geometry within the Geant4 world. The tool then provides additional useful
*        methods for calculating and returning information useful to various tracking
*        tasks.
*
* @author Tracy Usher
*
*/

#include "GlastSvc/Reco/IPropagator.h"
#include "ParticleTransporter.h"
#include "CLHEP/Matrix/Vector.h"
#include "G4Generator/IG4GenErrorSvc.h"
#include "G4Generator/G4GenException.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

class G4PropagationTool : public ParticleTransporter, public AlgTool, virtual public IPropagator
{
 public: 

    G4PropagationTool(const std::string& type, const std::string& name, const IInterface* parent);

    /// do stuff when set up by the ToolSvc
    virtual StatusCode initialize();

    virtual ~G4PropagationTool();

    //! Methods to intialize propagator before stepping 
    //! Tracking from intitial point and direction
    virtual void setStepStart(const Point& startPos, const Vector& startDir);
	//! Tracking from initial parameters
    virtual void setStepStart(const Event::TkrTrackParams& trackPar, double z, bool upwards);

    //! Takes a step of distance given by arcLen
    virtual void step(double arcLen);

    //! Number of steps taken
    virtual int getNumberSteps() const {return ParticleTransporter::getNumberSteps();}

    //! Methods to return information after the step
    //! Return position after stepping, arcLen can be less than step taken
    virtual Point getPosition(double arcLen = -1.) const;
    //! returns the arc length for this step
    virtual double getArcLen() const {return ParticleTransporter::getTotalArcLen();}
    //! return position at the center of the step in the volume
    virtual Point getStepPosition(int stepIdx = -1) const;
    //! returns the arc length for this step
    virtual double getStepArcLen(int stepIdx = -1) const;

    //! Return track parameters after stepping, arcLen can be less than step taken
    virtual Event::TkrTrackParams getTrackParams(double arcLen   = -1.,
                                                 double momentum = 1.,
                                                 bool   forward  = true)    const;
    //! Return multiple scattering matrix after stepping, arcLen can be less than step taken
    virtual HepMatrix getMscatCov(double arcLen   = -1.,
                                  double momentum =  1.,
                                  bool   forward  =  true) const;

    //! Return volume identifer after stepping
    virtual idents::VolumeIdentifier getVolumeId(double arcLen = -1.) const;
    //! Return volume identifer at given step index
    virtual idents::VolumeIdentifier getStepVolumeId(int stepIdx = -1) const;

    //! Return radiation lengths traversed
    virtual double getRadLength(double arcLen = -1.) const;
    //! Return radiation lengths for this step
    virtual double getStepRadLength(int stepIdx = -1.) const;

    //! Return the number of sensitive planes crossed to end of step
    virtual int    getNumSensePlanesCrossed() const;

    //! Is the point at the end of stepping inside an Active Area? 
    virtual double isInsideActArea()   const;
    virtual double isInsideActLocalX() const;
    virtual double isInsideActLocalY() const;

    //! Is the plane at the end of stepping an X or a Y plane?
    virtual bool   isXPlane()          const;
    virtual bool   isYPlane()          const;

    //! Is the strip at the end of stepping live?
    virtual bool   isStripLive()       const; 

    /// dump current status, to the stream
    virtual void printOn(std::ostream& str=std::cout )const;
  
    private:
    /// for use by isXPlane() and isYPlane()
    bool isTkrPlane(int& view) const;

    /// a pointer to the Error handling service
    IG4GenErrorSvc* m_ErrorSvc;

    /// Initial parameters
    Event::TkrTrackParams m_trackPar;
    double                m_zCoord;

    /// Maximum momentum to calculate ms elements
    double                m_maxMomentum;
};

//static ToolFactory<G4PropagationTool> g4prop_factory;
//const IToolFactory& G4PropagationToolFactory = g4prop_factory;
DECLARE_TOOL_FACTORY(G4PropagationTool);

G4PropagationTool::G4PropagationTool(const std::string& type, const std::string& name, const IInterface* parent) :
  ParticleTransporter(0,0),
  AlgTool(type, name, parent)
{
  // Purpose and Method:  Gaudi'ized constructor for the tool 
  // Inputs:  Stuff important to Gaudi
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  //Declare additional interface
  declareInterface<IPropagator>(this);

  // Make maximum momentum for calculation of ms errors a controllable parameter
  declareProperty("MaxMomentumForMS", m_maxMomentum = 50000.);

  //Initialize the track parameters
  m_trackPar = Event::TkrTrackParams();
  m_zCoord   = 0.;
}

StatusCode G4PropagationTool::initialize()
{
  // Purpose and Method:  Gaudi initialization routine. 
  // Inputs:  None
  // Outputs:  Gaudi StatusCode
  // Dependencies: None
  // Restrictions None 

  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "G4 Propagation Tool initializing now" << endreq;
  log << MSG::INFO << "G4Propagator tag v2r4p4 is initializing" << endreq;

  IService* iService = 0;
  if( serviceLocator()->getService( "G4GeometrySvc", iService, true).isFailure() ) 
  {
    log << MSG::ERROR << "Could not find the G4 Geometry Service!" << endreq;
    return StatusCode::FAILURE;
  }

  IG4GeometrySvc* gsv = dynamic_cast<IG4GeometrySvc*>(iService);

  // Error handling service 
  iService = 0;
  if (service("G4GenErrorSvc",iService, true).isFailure())
  {
      log << MSG::ERROR << "Could not find G4GenErrorSvc"<<endreq ;
      return StatusCode::FAILURE ;
  }

  m_ErrorSvc = dynamic_cast<IG4GenErrorSvc*>(iService);

  initExceptionHandler();

  setGeometrySvc(gsv);

  setTransportationManager(gsv->getTransportationManager());

  log << MSG::INFO << "G4 Propagation Tool ready" << endreq;

  return StatusCode::SUCCESS;
}

G4PropagationTool::~G4PropagationTool()
{
}

//! Tracking from intitial point and direction
void G4PropagationTool::setStepStart(const Point& startPos, const Vector& startDir)
{
    // Purpose and Method: Initializes to the starting point and direction. This
    //                     will serve to set the initial volume at the start point
    // Inputs: The starting point and direction 
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None
    m_trackPar(1) = startPos.x();
    m_trackPar(2) = startDir.x() / startDir.z();
    m_trackPar(3) = startPos.y();
    m_trackPar(4) = startDir.y() / startDir.z();

    m_trackPar(1,1) = 1.;
    m_trackPar(1,2) = 0.;
    m_trackPar(1,3) = 0.;
    m_trackPar(1,4) = 0.;
    m_trackPar(2,2) = 1.;
    m_trackPar(2,3) = 0.;
    m_trackPar(2,4) = 0.;
    m_trackPar(3,3) = 1.;
    m_trackPar(3,4) = 0.;
    m_trackPar(4,4) = 1.;

    m_zCoord   = startPos.z();

    setInitStep(startPos, startDir);

    return;
}

//! Tracking from initial parameters
void G4PropagationTool::setStepStart(const Event::TkrTrackParams& trackPar, double z, bool upwards)
{
    // Purpose and Method: Initializes to the starting point and direction. This
    //                     will serve to set the initial volume at the start point
    // Inputs: The starting point and direction 
    // Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats: None

    // Save parameters 
    m_trackPar = trackPar;
    m_zCoord   = z;

    // Create initial position and direction from said parameters
    double x_slope = trackPar.getxSlope();   
    double y_slope = trackPar.getySlope(); 
    Vector dir_ini = Vector(-x_slope, -y_slope, -1.).unit();
	if(upwards) dir_ini = Vector(x_slope, y_slope, 1.).unit();

    double x0      = trackPar.getxPosition();
    double y0      = trackPar.getyPosition();
    Point  x_ini(x0,y0,z);

    // Initialize the actual propagator...
    setInitStep(x_ini, dir_ini);

    return;
}

//! Takes a step of distance given by arcLen
void G4PropagationTool::step(double arcLen)
{
    try
    {
        transport(arcLen);
    }
    catch(G4GenException& e)
    {
        MsgStream log(msgSvc(), name());
        std::string exceptString = std::string("G4Exception - ") + std::string(e.what());
        
        log << MSG::WARNING << exceptString << std::endl; 

        throw std::runtime_error(exceptString);
    }
    catch(std::runtime_error& e)
    {
        MsgStream log(msgSvc(), name());

        log << MSG::WARNING << e.what() << endreq;

        throw e;
    }
    catch(...)
    {
        MsgStream log(msgSvc(), name());
        std::string exceptString("unknown exception caught");
        
        log << MSG::WARNING << exceptString << std::endl; 

        throw std::runtime_error(exceptString);
    }

    return;
}

//! Methods to return information after the step
//! Return position after stepping, arcLen can be less than step taken
Point G4PropagationTool::getPosition(double arcLen) const
{
    // Purpose and Method: Returns the position of the track after stepping arcLen distance
    //                     If arcLen is negative or more than stepped arcLen, returns position
    //                     at end of step.
    // Inputs:  None
    // Outputs:  a Point representing the final position
    // Dependencies: None
    // Restrictions and Caveats: None

    Point final(0.,0.,0.);

    if (getNumberSteps() > 0)
    {
        double totArcLen = getTotalArcLen();

        if (arcLen < 0. || arcLen >= totArcLen)
        {
            G4ThreeVector  stopPoint = getLastStep().GetEndPoint();

            final = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
        }
        else
        {
            G4ThreeVector  startPoint = getStepStart()->GetEntryPoint();
            G4ThreeVector  startDir   = getStepStart()->GetDirection();
            G4ThreeVector  stopPoint  = startPoint + arcLen * startDir;

            final = Point(stopPoint.x(),stopPoint.y(),stopPoint.z());
        }
    }

    return final;
}
Point G4PropagationTool::getStepPosition(int stepIdx) const
{
    // Purpose and Method: Returns the position of the track after stepping arcLen distance
    //                     If arcLen is negative or more than stepped arcLen, returns position
    //                     at end of step.
    // Inputs:  None
    // Outputs:  a Point representing the final position
    // Dependencies: None
    // Restrictions and Caveats: None

    int idx = stepIdx;

    if (idx < 0 || idx >= getNumberSteps()) idx = getNumberSteps() - 1;

    G4ThreeVector stepPoint = getStep(idx).GetEndPoint();

    return Point(stepPoint.x(),stepPoint.y(),stepPoint.z());
}

//! Return track parameters after stepping, arcLen can be less than step taken
Event::TkrTrackParams G4PropagationTool::getTrackParams(double arcLen,
                                                        double momentum,
                                                        bool   forward) const
{
    // Purpose and Method: Returns the track parameters after stepping arcLen distance
    //                     If arcLen is negative or more than stepped arcLen, returns position
    //                     at end of step.
    // Inputs:  None
    // Outputs:  a Point representing the final position
    // Dependencies: None
    // Restrictions and Caveats: None

    // The eventual result
    Event::TkrTrackParams trackPar;

    // Step taken along z axis
    double relDeltaZ = arcLen * getStartDir().z();

    // Look for exceeding limits
    if (arcLen < 0. || arcLen > getTotalArcLen())
    {
        relDeltaZ = getTotalArcLen() * getStartDir().z();
    }

    // Build a transport matrix
    HepMatrix F(4,4,1);
    F(1,2) = relDeltaZ;
    F(3,4) = relDeltaZ;

    // Vector representing the track parameters
    HepVector TrkP(4);
    TrkP(1) = m_trackPar(1);
    TrkP(2) = m_trackPar(2);
    TrkP(3) = m_trackPar(3);
    TrkP(4) = m_trackPar(4);

    TrkP = F * TrkP;

    trackPar(1) = TrkP(1);
    trackPar(2) = TrkP(2);
    trackPar(3) = TrkP(3);
    trackPar(4) = TrkP(4);

    // Matrix representing current covariance matrix
    HepMatrix cov(4,4,0);
    cov(1,1) = m_trackPar(1,1);
    cov(1,2) = m_trackPar(1,2);
    cov(1,3) = m_trackPar(1,3);
    cov(1,4) = m_trackPar(1,4);
    cov(2,1) = m_trackPar(2,1);
    cov(2,2) = m_trackPar(2,2);
    cov(2,3) = m_trackPar(2,3);
    cov(2,4) = m_trackPar(2,4);
    cov(3,1) = m_trackPar(3,1);
    cov(3,2) = m_trackPar(3,2);
    cov(3,3) = m_trackPar(3,3);
    cov(3,4) = m_trackPar(3,4);
    cov(4,1) = m_trackPar(4,1);
    cov(4,2) = m_trackPar(4,2);
    cov(4,3) = m_trackPar(4,3);
    cov(4,4) = m_trackPar(4,4);

    HepMatrix Qmaterial = getMscatCov(arcLen, momentum, forward);

    if (forward) cov = (F * (cov * F.T())) + Qmaterial;
    else         cov = (F * (cov + Qmaterial) * F.T());

    trackPar(1,1) = cov(1,1);
    trackPar(1,2) = cov(1,2);
    trackPar(1,3) = cov(1,3);
    trackPar(1,4) = cov(1,4);
    trackPar(2,1) = cov(2,1);
    trackPar(2,2) = cov(2,2);
    trackPar(2,3) = cov(2,3);
    trackPar(2,4) = cov(2,4);
    trackPar(3,1) = cov(3,1);
    trackPar(3,2) = cov(3,2);
    trackPar(3,3) = cov(3,3);
    trackPar(3,4) = cov(3,4);
    trackPar(4,1) = cov(4,1);
    trackPar(4,2) = cov(4,2);
    trackPar(4,3) = cov(4,3);
    trackPar(4,4) = cov(4,4);

    return trackPar;
}

double G4PropagationTool::getStepArcLen(int stepIdx) const
{
  // Purpose and Method:  Return the step length in the volume at stepIdx 
  // Inputs:  the index into the list of volumes traversed
  // Outputs:  the arc length in the selected volume
  // Dependencies: None
  // Restrictions None 

    int idx = stepIdx;

    if (idx < 0 || idx >= getNumberSteps()) idx = getNumberSteps();

    return getStep(stepIdx).GetArcLen();
}

//! Return track covariance matrix after stepping, arcLen can be less than step taken
//! *** NOTE *** this only has meaning if initial matrix given at setStepStart
/*
Event::TkrFitMatrix G4PropagationTool::getTrackCov(int propDir, double momentum, double arcLen) const
{
    // Purpose and Method: Returns the track parameters after stepping arcLen distance
    //                     If arcLen is negative or more than stepped arcLen, returns position
    //                     at end of step.
    // Inputs:  None
    // Outputs:  a Point representing the final position
    // Dependencies: None
    // Restrictions and Caveats: None

    double relDeltaZ = arcLen * getStartDir().z();

    if (arcLen < 0. || arcLen > getLastStep().GetArcLen())
    {
        relDeltaZ = getLastStep().GetArcLen() * getStartDir().z();
    }

    Event::TkrFitMatrix F(relDeltaZ);

    Event::TkrFitMatrix Qmaterial = getMscatCov(momentum, arcLen);
    Event::TkrFitMatrix cov;

    if (propDir < 0) cov = (F * (cov * F.T())) + Qmaterial;
    else             cov = (F * (cov + Qmaterial) * F.T());

    return cov;
}
*/
//! Return multiple scattering matrix after stepping, arcLen can be less than step taken
HepMatrix G4PropagationTool::getMscatCov(double arcLenIn, double momentum, bool) const
{
    // Purpose and Method: Calculates Q, the 4x4 covariance matrix due to multiple
    //                     scattering
    // Inputs: The total momentum of the track and the arc length over which to
    //         calculate
    // Outputs:  a 4x4 HepMatrix  
    // Dependencies: To get any answer, must have already stepped
    // Restrictions and Caveats: None

    double scat_dist  = 0.;
    double scat_angle = 0.; 
    double scat_covr  = 0.; 
    double dist       = 0.;
    double p33        = 0.;
    double p34        = 0.;
    double p44        = 0.;
    double arcLen     = arcLenIn;

    // Really nothing to do unless under the incoming momentum limit
    if (momentum < m_maxMomentum)
    {
        double x0sTotal   = 0.;    /// New stuff!

        if (arcLen < 0) arcLen = getTotalArcLen();

        ConstStepPtr stepPtr = getStepStart();

        while(stepPtr < getStepEnd())
        {
            TransportStepInfo  curStep = *stepPtr++;

            const G4VPhysicalVolume* pCurVolume = curStep.GetVolume();
            G4Material* pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();

            double radLengths = pMaterial->GetRadlen();
            double x0s        = 10000000.;
            double s_dist     = curStep.GetArcLen();
            double s_distp    = s_dist;

            if (radLengths > 0.) x0s = s_dist / radLengths;

            if(dist+s_dist > arcLen) { // pro-rate the last step: s_distp
                if(s_dist > 0) x0s *= (arcLen - dist)/s_dist;
                s_distp = arcLen - dist; 
            }

            if(x0s != 0) {
                ///double ms_Angle = 14.0*sqrt(x0s)*(1+0.038*log(x0s))/momentum; //MeV
                double ms_Angle = 13.6*sqrt(x0s)/momentum; //MeV  /// New Stuff!!
                double ms_Dst  = (arcLen - dist - s_distp)*ms_Angle; // Disp. over remaining traj
                double ms_sDst = s_distp*ms_Angle/1.7320508; // Disp. within step
                double ms_Dist = ms_Dst*ms_Dst + ms_sDst*ms_sDst;

                x0sTotal += x0s;   /// New Stuff!

                scat_dist  += ms_Dist;
                scat_angle += ms_Angle*ms_Angle;
                scat_covr  += sqrt(ms_Dist)*ms_Angle;		  
            }
            dist += s_dist;
            if(dist >= arcLen ) break;
        }

        /// New Stuff!!!!!
        if (x0sTotal > 0.)
        {
            static const double oneOverSqrt3 = 1.; // / sqrt(3.);
            // This returns the rms scattering angle in the plane of scatter
            double ms_Angle = 13.6*sqrt(x0sTotal)*(1+0.038*log(x0sTotal))/momentum; //MeV  //// New Stuff!!

            // The 1/sqrt(3) averages out the deviations over the step length
            scat_angle = oneOverSqrt3*oneOverSqrt3*ms_Angle*ms_Angle;

            // Try this to see if we get consistency
            double ms_dist = oneOverSqrt3 * arcLen * ms_Angle;
            scat_dist  = ms_dist*ms_dist;
            scat_covr  = ms_dist*ms_Angle;
        }
        else
        {
            scat_angle = 0.;
        }

        Vector startDir = getStartDir();
        double slopeX    = startDir.x()/startDir.z(); 
        double slopeY    = startDir.y()/startDir.z();
        // "norm_term" is actually 1 / cos(theta)^2 
        double norm_term = 1. + slopeX*slopeX + slopeY*slopeY;

        // The below taken from KalParticle (by Bill Atwood) in order to match results
        // Calc, the matrix elements (see Data Analysis Tech. for HEP, Fruhwirth et al)
        p33 = (1.+slopeX*slopeX)*norm_term;
        p34 = slopeX*slopeY*norm_term;
        p44 = (1.+slopeY*slopeY)*norm_term; 

        //Go from arc-length to Z - multiply by cos(theta) (or squared as case may be)
        scat_dist /= norm_term;
        scat_covr /= 2. * sqrt(norm_term);   // I do not know why I need the factor of two here...
    }

    HepMatrix cov(4,4,0);
    cov(1,1) = scat_dist*p33;
    cov(2,2) = scat_angle*p33; 
    cov(3,3) = scat_dist*p44;
    cov(4,4) = scat_angle*p44;
    cov(1,2) = cov(2,1) = -scat_covr*p33;
    cov(1,3) = cov(3,1) =  scat_dist*p34;
    cov(1,4) = cov(2,3) = cov(3,2) = cov(4,1) = -scat_covr*p34;
    cov(2,4) = cov(4,2) =  scat_angle*p34;
    cov(3,4) = cov(4,3) = -scat_covr*p44; 
  
    return cov;
}

//! Return radlen for this step
double G4PropagationTool::getStepRadLength(int stepIdx) const
{
    // Purpose and Method: Returns the arc length subtended at the end of tracking
    // Inputs:  None
    // Outputs:  a double representing the total arc length 
    // Dependencies: None
    // Restrictions and Caveats: None

    int idx = stepIdx;

    if (idx < 0 || idx >= getNumberSteps()) idx = getNumberSteps();

    const G4VPhysicalVolume* pCurVolume = getStep(idx).GetVolume();
    G4Material*        pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();

    double matRadLen = pMaterial->GetRadlen();
    double x0s       = 0.;
    double s_dist    = getStep(idx).GetArcLen();

    if (matRadLen > 0.) x0s = s_dist / matRadLen;

    return x0s;
}


//! Return radiation lengths traversed
double G4PropagationTool::getRadLength(double arcLenIn) const
{
    // Purpose and Method: Returns the arc length subtended at the end of tracking
    // Inputs:  None
    // Outputs:  a double representing the total arc length 
    // Dependencies: None
    // Restrictions and Caveats: None

    double radLen  = 0.;
    double totDist = 0.;
    double arcLen  = arcLenIn;

    if (arcLen < 0) arcLen = getTotalArcLen();

    //Set up iterator for stepping through all the layers
    ConstStepPtr stepPtr = getStepStart();

    while(stepPtr < getStepEnd())
    {
        TransportStepInfo  curStep = *stepPtr++;

        const G4VPhysicalVolume* pCurVolume = curStep.GetVolume();
        G4Material*              pMaterial  = pCurVolume->GetLogicalVolume()->GetMaterial();


        double matRadLen = pMaterial->GetRadlen();
        double x0s       = 0.;
        double s_dist    = curStep.GetArcLen();

        if (matRadLen > 0.) x0s = s_dist / matRadLen;

        // Check that next step is within max arc length
        if(totDist+s_dist > arcLen) 
        { 
            // pro-rate the last step: s_distp
            if(s_dist > 0) x0s *= (arcLen - totDist)/s_dist;
            s_dist = arcLen - totDist;
        }

        if (matRadLen > 0.) radLen += x0s;

        if ((totDist += s_dist) > arcLen) break;
    }

    return radLen;
}

//! Return the number of sensitive planes crossed to end of step
int G4PropagationTool::getNumSensePlanesCrossed() const
{
    // Purpose and Method: Returns the number of planes crossed. Tricky, first
    // volume is the //start point. The first swim is to this boundary, so the
    // volume is repeated. So, number of boundaries is size of vector - 2
    // Inputs:  None
    // Outputs:  Integer number of planes crossed
    // Dependencies: None
    // Restrictions and Caveats: None

    return getNumberSteps() > 1 ? getNumberSteps() - 1 : 0;
}

//! Is the point at the end of stepping inside an Active Area? 
double G4PropagationTool::isInsideActArea()   const
{
    // Purpose and Method:  Returns the distance to the nearest boundary of the
    //                      nearest sensitive volume. 
    // Inputs:  None
    // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
    // Dependencies: None
    // Restrictions and Caveats: None

    //  G4VPhysicalVolume* pCurVolume = getLastStep()->GetVolume();

    double dist = insideActiveArea();

    return dist;
}

double G4PropagationTool::isInsideActLocalX() const
{
    // Purpose and Method:  Returns the distance to the nearest X boundary of the
    //                      nearest sensitive volume. 
    // Inputs:  None
    // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
    // Dependencies: None
    // Restrictions and Caveats: None

    return insideActiveLocalX();
}

double G4PropagationTool::isInsideActLocalY() const
{
    // Purpose and Method:  Returns the distance to the nearest X boundary of the
    //                      nearest sensitive volume. 
    // Inputs:  None
    // Outputs:  distance to the edge, >=0 if inside, < 0 if outside
    // Dependencies: None
    // Restrictions and Caveats: None

    return insideActiveLocalY();
}

// is the volume at the end of the step a tracker plane?
bool   G4PropagationTool::isTkrPlane(int& view)        const
{
    // Purpose and Method:  Determines if current stop point is a tracker plane
    // Inputs:  None
    // Outputs:  a bool, true if current point is a tracker plane, false otherwise 
    // Dependencies: None
    // Restrictions and Caveats: None

//**    G4VPhysicalVolume* pCurVolume = getVolume(getLastStep().GetCoords(), true);

    idents::VolumeIdentifier id = constructId(getLastStep().GetEndPoint(), getLastStep().GetDirection());

    // This test is specifically for the tracker.
    // Return "false" otherwise
    bool  inTower    = id[0]==0;
    bool  inTkr      = id[3]==1;
    int   botTop     = id[6];
    
    view = id[5];

    if( !inTower || !inTkr || id.size()<7 || botTop>1 ) return false;
    return true;
}

//! Is the plane at the end of stepping an X or a Y plane?
bool   G4PropagationTool::isXPlane()          const
{
    // Purpose and Method:  Determines if current stop point is an "X" plane
    // Inputs:  None
    // Outputs:  a bool, true if current point is an "X" plane, false otherwise 
    // Dependencies: None
    // Restrictions and Caveats: None

    int view;
    if (!isTkrPlane(view)) return false;

    return (view==0);
}

bool   G4PropagationTool::isYPlane()          const
{
    // Purpose and Method:  Determines if current stop point is an "Y" plane
    // Inputs:  None
    // Outputs:  a bool, true if current point is an "Y" plane, false otherwise 
    // Dependencies: None
    // Restrictions and Caveats: None

    int view;
    if (!isTkrPlane(view)) return false;

    return (view==1);
}

//! Is the strip at the end of stepping live?
bool   G4PropagationTool::isStripLive()       const
{
    return true;
}


//! Return volume identifer after stepping
idents::VolumeIdentifier G4PropagationTool::getVolumeId(double arcLen) const
{
    // Purpose and Method:  Returns the volume identifier at end of final step
    // Inputs:  None
    // Outputs:  an idents::VolumeIdentifier 
    // Dependencies: None
    // Restrictions and Caveats: None

    ConstStepPtr stepPtr = getStepAtArcLen(arcLen);

//    G4VPhysicalVolume* pCurVolume = getVolume(stepPtr->GetCoords(), true);

    return constructId(stepPtr->GetEndPoint(), stepPtr->GetDirection(), true);
}

//! Return volume identifer after stepping
idents::VolumeIdentifier G4PropagationTool::getStepVolumeId(int stepIdx) const
{
    // Purpose and Method:  Returns the volume identifier at requested step
    // Inputs:  None
    // Outputs:  an idents::VolumeIdentifier 
    // Dependencies: None
    // Restrictions and Caveats: None

    int idx = stepIdx;

    if (idx < 0 || idx >= getNumberSteps()) idx = getNumberSteps();

//    G4VPhysicalVolume* pCurVolume = getVolume(getStep(stepIdx).GetCoords(), true);

    return constructId(getStep(stepIdx).GetEndPoint(), getStep(stepIdx).GetDirection(), true);
}


/// dump current status, to the stream
void G4PropagationTool::printOn(std::ostream& str)const
{
  // Purpose and Method:  Call through for printing tracking information 
  // Inputs:  a stream pointer
  // Outputs:  None
  // Dependencies: None
  // Restrictions None 

    str << "\n";

    printStepInfo(str);

    return;
}

//Header: /cvs/glastsim/flux/FluxSource.cxx,v 1.20 1998/12/18 23:42:02 hsg Exp $

#include "FluxSvc/FluxSource.h"

#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "xml/Dom.h"

#include "FluxSvc/SpectrumFactoryTable.h"

#include "FluxSvc/SimpleSpectrum.h"

#include "CLHEP/Random/RandFlat.h"

#include "geometry/CoordTransform.h"
#include "geometry/Ray.h"
#include "geometry/Box.h"

#include "Orbit.h"

#include "FluxException.h" // for FATAL_MACRO

double  FluxSource::s_backoff;
double  FluxSource::s_radius=1.0;

// display/command interfaces...
#include <algorithm>
#include <list>
using std::min;
using std::max;


FluxSource::FluxSource(Spectrum* aSpec, double aFlux)
: EventSource(aFlux),  m_spectrum(0),
m_maxEnergy(100.),  // note defualt maximum kinetic energy
_minCos(-0.4f), _maxCos(1.0f), _minPhi(0.0f), _maxPhi(2*M_PI),
m_rmin(0), m_rmax(1), _phi(0.0f), _theta(0.0f), m_pointtype(NOPOINT),
m_launch(NONE),m_frametype(EARTH), illumBox(0)
{
    s_backoff = 0.;
    spectrum(aSpec);
    useSpectrumDirection(); // and will use its direction generation
    setAcceptance();
  //  transformDirection();
}

FluxSource::FluxSource(const DOM_Element& xelem )
: EventSource (xelem), m_spectrum(0),
m_maxEnergy(100.),  // note defualt maximum kinetic energy
_minCos(-0.4f), _maxCos(1.0f), _minPhi(0.0f), _maxPhi(2*M_PI),
m_rmin(0), m_rmax(1), _phi(0.0f), _theta(0.0f), m_pointtype(NOPOINT), m_launch(NONE),m_frametype(EARTH), illumBox(0)
{
    static double d2r = M_PI/180.;
    
    Spectrum*   s = 0;
    std::string class_name;
    std::string source_params; 
    
    DOM_Element   spec = xml::Dom::findFirstChildByName(xelem, "spectrum");
    
    
    if (spec == DOM_Element()) {
        
        // source has no imbedded spectrum element: expect a name
        class_name = xelem.getAttribute("name").transcode();
        useSpectrumDirection(); // and will use its direction generation
        
    } else {
        
        // process spectrum element
        DOM_NodeList children = spec.getChildNodes();
        
        // First child element is type of spectrum
        DOM_Node    childNode = children.item(0);
        DOM_Element specType;
        
        if (childNode.getNodeType() == DOM_Node::ELEMENT_NODE) {
            specType = (DOM_Element &) childNode;
        }
        else specType = xml::Dom::getSiblingElement(childNode);
        
        DOMString   typeTagName = specType.getTagName();
        std::string spectrum_name = spec.getAttribute("name").transcode();
        std::string spectrum_frametype = spec.getAttribute("frame").transcode();

        if(spectrum_frametype == "galactic"){ m_frametype=GALAXY;
        }else if(spectrum_frametype == "glast"){ m_frametype=GLAST;
        }else if(spectrum_frametype == "earth"){ m_frametype=EARTH;
        }

        if (typeTagName.equals("particle")) s = new SimpleSpectrum(specType);
        else if (typeTagName.equals("SpectrumClass")) {
            // attribute "name" is the class name
            class_name = specType.getAttribute("name").transcode();
            source_params= specType.getAttribute("params").transcode();
        }
        else {
            // no, the tag itself
            class_name = typeTagName.transcode();
        }
        // second child element is angle
        DOM_Element angles = xml::Dom::getSiblingElement(specType);
        
        DOMString anglesTag = angles.getTagName();
        
        if (anglesTag.equals("solid_angle") ) {
            setLaunch(atof(angles.getAttribute("theta").transcode()) * d2r, 
                atof(angles.getAttribute("phi").transcode()) *d2r);
            setCosThetaRange(atof(angles.getAttribute("mincos").transcode()),
                atof(angles.getAttribute("maxcos").transcode())  );
        }
        else if (anglesTag.equals("direction") ) {
            setLaunch(atof(angles.getAttribute("theta").transcode()) * d2r, 
                atof(angles.getAttribute("phi").transcode()) *d2r);
        }
        else if (anglesTag.equals("use_spectrum") ) {
            useSpectrumDirection();
        }
        else if(anglesTag.equals("galactic_dir")){
            m_galb=atof(angles.getAttribute("b").transcode());
            m_gall=atof(angles.getAttribute("l").transcode());
            getGalacticDir(atof(angles.getAttribute("l").transcode()),
                atof(angles.getAttribute("b").transcode()) );
            m_launch=GALACTIC;
        }
        else {
            FATAL_MACRO("Unknown angle specification in Flux::Flux \""
                << anglesTag.transcode() << "\"" );
        }
        
        // third child element is optional launch spec
        DOM_Element launch = xml::Dom::getSiblingElement(angles);
        if(launch !=DOM_Element()) {
            
            DOMString launchTag = launch.getTagName();
            
            if(launchTag.equals("launch_point")){
                m_pointtype=SINGLE;
                LaunchType templaunch=m_launch;
                // assume that the launch direction was set by a direction!
                setLaunch(/*launchDir()*/m_launchDir, 
                    Point(atof(launch.getAttribute("x").transcode()),
                    atof(launch.getAttribute("y").transcode()),
                    atof(launch.getAttribute("z").transcode()) ) );
                m_launch=templaunch;
                
            }else if(launchTag.equals("patch")){
                
                m_pointtype = PATCH;
                LaunchType templaunch=m_launch;
                
                //setLaunch(launchDir().theta(), launchDir().phi(),
                
                setLaunch( /*m_launchDir*/(-launchDir()).theta(), /*m_launchDir*/(-launchDir()).phi(),
                    atof(launch.getAttribute("xmax").transcode()),
                    atof(launch.getAttribute("xmin").transcode()),
                    atof(launch.getAttribute("ymax").transcode()),
                    atof(launch.getAttribute("ymin").transcode()), 
                    atof(launch.getAttribute("zmax").transcode()),
                    atof(launch.getAttribute("zmin").transcode()) 
                    );
                
                m_launch=templaunch;
                
                
            }else {
                FATAL_MACRO("Unknown launch specification in Flux::Flux \""
                    << launchTag.transcode() << "\"" );
            }
        }
        
        
    }
    if( s==0) {
        //		std::vector<float> paramvec; parseParamList(source_params, paramvec);
        s = SpectrumFactoryTable::instance()->instantiate(class_name, source_params);
        
        if(s==0){
            FATAL_MACRO("Unknown Spectrum: "<< class_name);
            std::cerr << "List of known Spectrum classes:\n" ;
            std::list<std::string>list= SpectrumFactoryTable::instance()->spectrumList();
            for( std::list<std::string>::iterator i = list.begin(); i!=list.end(); ++i)
                std::cerr << "\t" << *i << std::endl;
            
            return;
        }
    }
    
    // finally set the spectrum object, and prepare to use it.
    FluxSource::spectrum(s);
    
    s_backoff = 0.;
    setAcceptance();

    //transformDirection();
    
}


FluxSource::~FluxSource()
{
    delete m_spectrum;
    if(illumBox) delete illumBox;
}

inline static double    sqr(double x) {return x*x; }


void FluxSource::setAcceptance()
{
    s_radius = sqrt(totalArea() / M_PI ) * 100;    // radius in cm
    if(!s_backoff) s_backoff = s_radius;
    
    switch (m_launch) {
        
    case NONE:  // cos theta range...
        EventSource::solidAngle( (_maxCos - _minCos)*(_maxPhi - _minPhi) );
        break;
    case DIRECTION: // single direction solid angle = 1.
    case POINT:         // treat as a direction - solid angle = 1.
    case GALACTIC:  //single direction, solid angle = 1.
    case SPECTRUM:
        // make sure that the rate calculation uses the spectrum's count
        EventSource::solidAngle( m_spectrum==0? 1.0 : m_spectrum->solidAngle());
        break;
    case SURFACE:
        EventSource::solidAngle( (_maxCos - _minCos) * (_maxPhi - _minPhi) );
        if (_maxCos == _minCos) EventSource::solidAngle(1.);  // treat as direction solid angle = 1
        break;
    case PATCHFIXED:
        EventSource::solidAngle(1.0); // single direction, solid angle = 1.
        break;
    }
}

void FluxSource::spectrum(Spectrum* s, double emax)
{
    if (emax > 0) {
        setMaxEnergy(emax);
        //m_rmax =  s->fraction(emax);
        std::cerr << "exercising obsolete function fraction" << std::endl;
    }
    
    m_spectrum = s;
    //const char* name = s->particleName();
}

FluxSource* FluxSource::event(double) 
{
    computeLaunch();
    return this;
    // could be a call-back
}

void FluxSource::randomLaunchPoint()
{
    HepRotation r_pln;
    
    double ly = m_launchDir.y(), lx = m_launchDir.x();
    if( lx !=0 || ly !=0 )r_pln.rotate(acos(m_launchDir.z()), 
        Vector(-ly, lx, 0.));
    
    // pick a random position on the planar section of a sphere through 
    // its midpoint
    double azimuth = RandFlat::shoot( 2*M_PI );
    double rad = s_radius*(sqrt(RandFlat::shoot()));
    
    // create two vectors to describe the particle launch: one to describe
    // the point in the plane perpendicular to the launch direction (within
    // the cross-section of the sphere containing the instrument) and a 
    //second to describe the distance along the normal between the launch 
    // point and that plane.
    Vector posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);
    
    // rotate these vectors according to the theta, phi specs
    CoordTransform  t(r_pln, Vector(0,0,0) ); // Vector(_box->center())
    t.transformVector(posLaunch);
    
    // define actual launch point
    m_launchPoint = (Point&)posLaunch - m_launchDir*s_backoff;
}



double FluxSource::flux(double time) const
{
    return enabled() ? std::max(m_spectrum->flux(time),0./* EventSource::flux()*/) : 0;
}

#if 1
double FluxSource::solidAngle() const
{
    //  return std::max(m_spectrum->solidAngle(), EventSource::solidAngle());
    return  EventSource::solidAngle();
}
#endif

void FluxSource::computeLaunch ()
{
    // set energy using the Spectrum object (scales momentum)
    // Note: since PEGS files crap out at some energ, the max energy must
    // be limited
    const double fudge=1.001; // in ncase max is only energy, round-off error
    double kinetic_energy;
    do {
        // kinetic_energy= (*spectrum())(RandFlat::shoot(m_rmin, m_rmax));
        //FIXME: make this a class variable
        kinetic_energy = spectrum()->energySrc( HepRandom::getTheEngine() );
    }    while (kinetic_energy > m_maxEnergy* fudge);
    
    // get the launch point and direction, according to the various strategies
    
    switch (m_launch) {
    case SPECTRUM:
        {
            // special option that gets direction from the spectrum object
            // note extra - sign since direction corresponds to *from*, not *to*
            std::pair<float,float> direction = spectrum()->dir(kinetic_energy,HepRandom::getTheEngine());
            double costh = direction.first;
            double sinth = sqrt(1.-costh*costh);
            double phi = direction.second;
            Vector v(cos(phi)*sinth, sin(phi)*sinth, costh);
            m_launchDir = -v;
            
            // extra rotation in case GLAST is not zenith pointing (beware, 
            // might be confusing)
            // keep x-axis perpendicular to zenith direction
            if (_theta != 0.0) m_launchDir.rotateX(_theta).rotateZ(_phi);
            
            //if(m_pointtype==NOPOINT){
            //randomLaunchPoint();
            //}
            break;
        }
    case GALACTIC:
        {
            //getGalacticDir should already have been called in the constructor
            
            getGalacticDir(m_gall, m_galb);
            break;
        }
    case SURFACE:
        {
            getSurfacePosDir();
            break;
        }
    case PATCHFIXED: {
        // Check for normal incidence...don't want to try to rotate about a 
        // zero vector
        if (_theta == 0.0) {
            // Just choose a random position within the rectangle
            double xInc = patchXmin + patchWidX * RandFlat::shoot();
            double yInc = patchYmin + patchWidY * RandFlat::shoot();
            double zInc = patchTop;
            m_launchPoint = Point(xInc, yInc, zInc);
        } else {
            double dist = FLT_MAX;
            const double distToSearch = 500.;
            do {
                randomLaunchPoint();   // pick a point on the sphere as usual
                // Test to see if this position pierces our patch
                Ray trialRay(m_launchPoint, m_launchDir);
                dist = illumBox->distanceToEnter(trialRay, distToSearch);
                // set the launchPoint to begin on our surface - the patch
                if (dist < FLT_MAX) m_launchPoint = trialRay.position(dist);
                
                // search til we find a position that satisfies our patch
            } while (dist >= FLT_MAX); 
        }
        break;
                     }
    case NONE:
        {
            double  costh = -RandFlat::shoot(_minCos, _maxCos);
            double  sinth = sqrt(1.-costh*costh);
            double  phi = RandFlat::shoot(_minPhi, _maxPhi);
            
            // extra rotation in case not zenith pointing (beware, might be
            // confusing)
            m_launchDir = Vector(cos(phi)*sinth, sin(phi)*sinth, costh);
            
            // keep x-axis perpendicular to zenith direction
            if (_theta != 0.0) m_launchDir.rotateX(_theta).rotateZ(_phi);
        }
        // fall through to...
    case DIRECTION:
        {
            //unused
        }
        // fall through to...
    case POINT:
        // here now that point and direction are set
        break;
        
    } // switch m_launch
    
    m_energy = kinetic_energy;
    
    
    
    switch (m_pointtype) {
    case NOPOINT:
        {
            randomLaunchPoint();
            break;
        }
        
    case SINGLE:
        {
            break;
        }
    case PATCH:
        {
            // Check for normal incidence...don't want to try to rotate about a 
            // zero vector
            if (_theta == 0.0) {
                // Just choose a random position within the rectangle
                double xInc = patchXmin + patchWidX * RandFlat::shoot();
                double yInc = patchYmin + patchWidY * RandFlat::shoot();
                double zInc = patchTop;
                m_launchPoint = Point(xInc, yInc, zInc);
            } else {
                double dist = FLT_MAX;
                const double distToSearch = 500.;
                do {
                    randomLaunchPoint();   // pick a point on the sphere as usual
                    // Test to see if this position pierces our patch
                    Ray trialRay(m_launchPoint, m_launchDir);
                    dist = illumBox->distanceToEnter(trialRay, distToSearch);
                    // set the launchPoint to begin on our surface - the patch
                    if (dist < FLT_MAX) m_launchPoint = trialRay.position(dist);
                    
                    // search til we find a position that satisfies our patch
                } while (dist >= FLT_MAX); 
            }
            break;
        }
    }
//   transformDirection(); 
    
}

std::string FluxSource::fullTitle () const
{
    return title();
}

std::string FluxSource::displayTitle () const
{
    std::strstream s;
    s << EventSource::displayTitle() << '(' << m_spectrum->title() ;
    s << ')' << '\0';
    std::string t(s.str());
    s.freeze(false);
    return t;
}

int FluxSource::eventNumber()const
{
    return 0;
}

void  FluxSource::setCosThetaRange(double minc, double maxc)
{
    // require _maxCos > _minCos for solid angle calculation
    _minCos   = min(minc,maxc);
    _maxCos   = max(minc,maxc);
    if(_minCos==_maxCos) {
        if(_minCos!=-1) _minCos-=0.001; else _maxCos +=0.001;
    }
    m_launch  = NONE;
    setAcceptance();
}

void  FluxSource::setPhiRange(double min_phi, double max_phi)
{
    // Assume angles are given in degrees
    _minPhi  = min_phi;
    _maxPhi  = max_phi;
    m_launch = NONE;
    setAcceptance();
}
void FluxSource::setLaunch(const Vector& dir, const Point& pos)
{
    //std::cout << "setting launch" << std::endl;
    m_launchDir   = dir;
    m_launchPoint = pos;
    m_launch     = POINT;
}
void FluxSource::useSpectrumDirection()
{
    m_launch = SPECTRUM;
}


void FluxSource::setLaunch(double theta, double phi)
{
    //std::cout << "setting launch" << std::endl;
    _theta = theta;
    _phi   = phi;
    Vector dir(sin(_theta)*cos(_phi),sin(_theta)*sin(_phi),cos(_theta));
    setLaunch(-dir); // minus due to z axis pointing UP!
}

void FluxSource::setLaunch(const Vector& dir)
{
    //std::cout << "setting launch" << std::endl;
    m_launchDir    = dir.unit();
    m_launch      = DIRECTION;
    _theta = (-dir).theta();
    _phi   = (-dir).phi();
    setAcceptance();
}

void FluxSource::setLaunch(double theta, double phi, double xMax, 
                           double xMin, double yMax, double yMin, 
                           double zTop, double zBot) {
    
    //std::cout << "setting launch" << std::endl;
    // fixed angle patch    
    // Just set up the area of illumination
    double epsilon = 1e-5;
    
    if (zTop < zBot) {
        WARNING_MACRO("zTop < zBot -- SWAPPING");
        std::swap(zTop, zBot);
    }
    
    patchTop = zTop;
    patchBottom = zBot;
    patchHeight = fabs(zTop - zBot);
    if(patchHeight < epsilon) patchHeight = 0.0;
    
    if (xMax < xMin) {
        WARNING_MACRO("patchXmax < patchXmin in Flux -- SWAPPING");
        std::swap(xMax, xMin);
    }
    
    if( yMax < yMin) {
        WARNING_MACRO("patchYmax < patchYmin in Flux -- SWAPPING");
        std::swap(yMax, yMin);
    }
    
    patchXmax = xMax;
    patchXmin = xMin;
    patchYmax = yMax;
    patchYmin = yMin;
    
    patchWidX = fabs(patchXmax - patchXmin);
    patchWidY = fabs(patchYmax - patchYmin);
    if(patchWidX < epsilon) patchWidX = 0.0;
    if(patchWidY < epsilon) patchWidY = 0.0;
    
    if( (patchHeight == 0.0) && ( (patchWidX == 0.0) || (patchWidY == 0.0) ) ) {
        FATAL_MACRO("zero area illumination\n");
    }
    if( (patchWidX == 0.0) && (patchWidY == 0) ) {
        FATAL_MACRO("zero area illumination!\n");
    }
    
    // Check if it's on-axis but the height of our box is non-zero
    if ((theta == 0.0) && (patchHeight != 0.0)) {
        WARNING_MACRO("on-axis, but zTop != zBot, setting zBot = zTop\n");
        patchHeight = 0.0;
        patchBottom = patchTop;
    }
    
    // If one of the dimensions is missing...create a small one
    // We want to use a box - so this forces the issue
    if (patchHeight <= 0.0) patchHeight = 0.0001;
    if (patchWidX <= 0.0) patchWidX = 0.0001;
    if (patchWidY <= 0.0) patchWidY = 0.0001;
    
    // Setup our box...used to test when we have chosen on position on 
    // the sphere that will work
    illumBox = new Box(patchWidX, patchWidY, patchHeight);
    Vector center((xMax + xMin)/2., (yMax+yMin)/2., (zBot+zTop)/2.);
    illumBox->move(center);
    
    // setup the fixed angles
    _theta = theta;
    _phi   = phi;
    Vector dir(sin(_theta)*cos(_phi),sin(_theta)*sin(_phi),cos(_theta));
    // minus dir since z axis is pointing up, see setLaunch(theta, phi)
    m_launchDir    = (-dir).unit();
    m_launch      = PATCHFIXED;
    
    setAcceptance();
}

void FluxSource::setLaunch(double xMax, double xMin, double yMax, double yMin,
                           double zTop, double zBot, bool fan)
{
    
    
    // Inputs: bounds of the incident box: xMax, xMin, yMax, yMin, zTop, zBot
    //         fan denotes whether or not this is a fan beam
    // setup illumination of some portion of the instrument - could be one side
    // or a box.  If illuminating one side - it will always be the +X side 
    // of GLAST
    
    double epsilon = 1e-5;
    
    sidePatch = false;
    fanBeam = false;
    double thMax = acos(_minCos);  // theta = [0, PI]
    double thMin = acos(_maxCos);
    
    if (zTop < zBot) {
        WARNING_MACRO("zTop < zBot -- SWAPPING");
        std::swap(zTop, zBot);
    }
    
    patchTop = zTop;
    patchBottom = zBot;
    patchHeight = fabs(zTop - zBot);
    if(patchHeight < epsilon) patchHeight = 0.0;
    
    if (xMax < xMin) {
        WARNING_MACRO("patchXmax < patchXmin in Flux -- SWAPPING");
        std::swap(xMax, xMin);
    }
    
    if( yMax < yMin) {
        WARNING_MACRO("patchYmax < patchYmin in Flux -- SWAPPING");
        std::swap(yMax, yMin);
    }
    
    patchXmax = xMax;
    patchXmin = xMin;
    patchYmax = yMax;
    patchYmin = yMin;
    
    patchWidX = fabs(patchXmax - patchXmin);
    patchWidY = fabs(patchYmax - patchYmin);
    if(patchWidX < epsilon) patchWidX = 0.0;
    if(patchWidY < epsilon) patchWidY = 0.0;
    
    // area of the sides of GLAST
    double aTop, aBot, aSide0, aSide1, aSide2, aSide3; 
    
    // Define the areas of the sides that may be hit by the beam.
    if( (patchWidX == 0.0) || (patchWidY == 0.0) ) {
        
        if( (patchHeight == 0.0) || ((patchWidX == 0.0) && (patchWidY == 0.0)) ) {
            FATAL_MACRO("zero area illumination\n");
        }
        
        // we're illuminating one side of the instrument
        sidePatch = true; 
        
        fanBeam = fan; // fan beam or not (aka polar beam)?
        
        if( patchWidY == 0.0) { 
            // always shoot into +x-axis, so swap the variables to reflect this
            double tempY = patchYmax;
            patchWidY = patchWidX;
            patchYmax = patchXmax;
            patchYmin = patchXmin;
            patchWidX = 0.0;
            patchXmin = tempY;
            patchXmax = tempY;
        }
        
        // compute the areas of the sides that may be illuminated
        aSide0 = patchHeight * patchWidY;  // +X axis side
        aSide1 = 0.0;
        aSide2 = 0.0;
        aSide3 = 0.0;
        aTop = 0.0;
        aBot = 0.0;
        
    } else { // illuminating more than just one side
        // Note this always refers to a polar beam, just not sidePatch
        aTop = patchWidX * patchWidY;
        aBot = patchWidX * patchWidY;
        aSide0 = patchHeight * patchWidY; // +X axis side
        aSide1 = patchHeight * patchWidX;
        aSide2 = patchHeight * patchWidY;
        aSide3 = patchHeight * patchWidX;
    }
    
    // Checking phi range for the case where we're illuminating one side of GLAST
    if( (sidePatch == true) && (fanBeam == false) ) { 
        //      WARNING_MACRO("illuminating one side with polar beam, PHI = (-PI, PI)");
        
    } else if ((sidePatch == true) && 
        ((_minPhi < (-M_PI/2.) * (1.+epsilon) ) ||
        (_maxPhi > (M_PI/2.) * (1.+epsilon)  )   )     )
    {  // fan Beam, side Patch
        // WARNING_MACRO("minPhi or maxPhi is out of range for illuminating one side with fan beam - reset to (-PI/2, PI/2)");
        _minPhi = -M_PI/2.;
        _maxPhi = M_PI/2.;
        
    } else if (sidePatch == false) { // illuminating more than one side of GLAST
        //      WARNING_MACRO("sidePatch == false:  require PHI = (-PI, PI)");
        _minPhi = -M_PI;
        _maxPhi = M_PI;
    }
    
    // total area to be illuminated
    double aSideTot = aSide0 + aSide1 + aSide2 + aSide3; 
    
    patchRange = (_minCos < 0 ? -1 : 1) * _minCos * _minCos -
        (_maxCos < 0 ? -1 : 1) * _maxCos * _maxCos;
    
    patchOffset = (_maxCos < 0 ? -1 : 1) * _maxCos * _maxCos;
    
    // cout << "patchRange, patchOffset = " << patchRange << " " << patchOffset << "\n";
    
    float wBot = 0.0, wTop = 0.0;
    
    if (sidePatch == true) {  
        // assumes we're hitting +x side, so no top or bottom
        wBot = 0.0;
        wTop = 0.0;
        
    } else {
        // Calculate geometric factor weights for Top and Bottom
        // based on integration of cos(theta) * d(Omega) about z-axis
        if ( (thMax <= M_PI / 2.) && (thMin < M_PI / 2.) ) {  
            // will only hit top, not bottom
            wBot = 0.0;
            wTop = _maxCos * _maxCos - _minCos * _minCos;
            // check for normal incidence
            if((thMax == 0) && (thMin == 0) ) wTop = 1.0;  
            
        } else if ( (thMax > M_PI / 2.) && (thMin < M_PI / 2.) ) {
            // could hit top or bottom
            wBot = _minCos * _minCos;
            wTop = _maxCos * _maxCos;
            
        } else { // (thMax > 90.) && (thMin >= 90.)  // only hit bottom, not top
            wBot = _minCos * _minCos - _maxCos * _maxCos;
            wTop = 0.0;
            if ( (thMax == M_PI) && (thMin == M_PI) ) wBot = 1.0; // bottom only
        }
    }
    
    // Calculate Weights for each side, based on integration
    // about z-axis of d(area) * d(Omega) =
    // (cos(phi) * sin(theta))*sin(theta) * d(phi) * d(theta)
    // range of phi is [-PI, PI] for sidePatch = false
    float wSide;
    
    if(sidePatch == true) {
        wSide = 1.0;
    } else {
        wSide = ( (thMax - thMin) - 0.5 * ( sin(2. * thMax) - sin(2. * thMin) ) ) 
            / M_PI;
    }
    
    Fratio = wSide * aSideTot / (wSide * aSideTot + wTop * aTop + wBot * aBot);
    // cout << "Fratio = " << Fratio << "\n";
    // cout << "wSide, wTop, wBot, aSideTot = " << wSide << " " << wTop << " " << wBot << " " << aSideTot << "\n";
    // cout << "aTop, aBot, aSide0, aSide1, aSide2, aSide3 = " << aTop << " " << aBot << " " << aSide0 << " " << aSide1 << " " << aSide2 << " " << aSide3 << "\n";
    // cout << "thMax, thMin, _minPhi, _maxPhi = " << thMax << " " << thMin << " " << _minPhi << " " << _maxPhi << "\n";
    
    m_launch = SURFACE;
    setAcceptance();
}

void FluxSource::getSurfacePosDir() {
    // "patch" uniform isotropic illumination schemes:
    //               more than one side:  4*PI, or delimited by polar cones
    //               one side (or part of one side)
    //               top / bottom (or part of top / bottom)
    //
    // Range of azimuthal angle, PHI:
    //   if more than one side is illuminated, range = (-PI,PI)
    //   if rectangle on one side (+X) is illuminated, range = (-PI/2,PI/2)
    //
    // One-side illumination:
    //   GLAST is quadrilaterally symmetric, so we assume the +X axis side.
    //   Side illumination has two options:  a fan beam or polar beam.
    //   In the fan-beam case, the polar range (w/ +Z-axis = symmetry axis)
    //   and azimuthal ranges can be delimited.  Polar-beam illumination on
    //   a side is realized by making the +X-axis the symmetry axis.
    //
    // 4*PI, delimited cones (> one side), or top/bottom illumination:
    //   Presently, the only option is with +Z-axis = symmetry axis
    //   with PHI range = (-PI,PI)
    //
    // Note that a rectangular-solid volume, or rectangular surface,
    //   interior to the GLAST instrument can be illuminated.
    double xInc, yInc, zInc;
    double cosTh, sinTh, phi;
    int sideNum;
    if (RandFlat::shoot() < Fratio) { // Hit a Side
        double cosE, sinE;
        
        if ( (sidePatch == true) && (fanBeam == false) ) {   
            // Polar Beam and Patch on +X side
            // Polar Beam, so theta = [ThetaMin, ThetaMax] about +X axis 
            // (not +Z axis)
            double randNum = RandFlat::shoot() * patchRange + patchOffset;
            cosE = ( (randNum < 0) ? -1 : 1) * sqrt( fabs (randNum) );
            sinE = sqrt( 1 - cosE * cosE);
            
            // Xi = [0, 2M_PI];  can be anything
            double Xi = RandFlat::shoot(0., 2*M_PI);  
            
            cosTh = sinE * sin(Xi);    // law of cosines
            sinTh = sqrt(1 - cosTh * cosTh);
            
            double anySidePhi = (Xi < 0 ? -1 : 1) * acos(cosE / sinTh);
            phi = anySidePhi;   // phi = [ -PI/2, PI/2]  always hitting +X side
            
            // already know that sidePatch == true, so we must hit +X side
            sideNum = 0;  
            
        } else {
            // either Fan Beam Patch on +X axis side OR illuminating more than
            // one side with a Polar Beam, with azimuthal coordinate unrestricted
            bool PhiDONE;
            do {
                PhiDONE = false;
                cosE = sqrt(RandFlat::shoot());  
                // get cosE = sqrt[0,1] which corresponds to [0, PI/2]
                // cosTheta isn't restricted until the conditional at the end 
                // of the do-while loop
                sinE = sqrt( 1 - cosE * cosE);
                
                double Xi = RandFlat::shoot(-M_PI, M_PI);
                cosTh = sinE * sin(Xi);
                sinTh = sqrt(1 - cosTh * cosTh);
                
                double rand = RandFlat::shoot(-1, 1);
                
                double anySidePhi = (rand < 0 ? -1 : 1) * acos(cosE / sinTh);
                
                sideNum = int(4 * RandFlat::shoot()); // choose which side to hit
                if(sidePatch == true) sideNum = 0;  // always hit +x side
                phi = anySidePhi + (M_PI / 2.) * sideNum;
                
                if(fanBeam == true) {  // is Phi within bounds?
                    if ( (phi >= _minPhi) && (phi <= _maxPhi) ) PhiDONE = true;
                } else {
                    PhiDONE = true;
                }
                
            } while( (cosTh > _maxCos) || (cosTh < _minCos) || (!PhiDONE) );
            
        }
        // choose incident position on side
        zInc = patchBottom + patchHeight * RandFlat::shoot();
        
        if ( (sideNum % 2) == 0) { // hit X side
            if(sideNum == 0) {
                xInc = patchXmax;
            } else {
                xInc = patchXmin;
            }
            yInc = patchYmin + patchWidY * RandFlat::shoot();
        } else {  // Hit Y side
            xInc = patchXmin + patchWidX * RandFlat::shoot();
            if (sideNum == 1) {
                yInc = patchYmax;
            } else {
                yInc = patchYmin;
            }
        }
    } else {
        
        // Hit top or bottom, within delimited polar cones; azimuthal 
        // coordinate unrestricted
        phi = RandFlat::shoot(_minPhi, _maxPhi);
        double ranN = RandFlat::shoot() * patchRange + patchOffset;
        cosTh = (ranN > 0 ? 1 : -1) * sqrt( fabs ( ranN ) );
        sinTh = sqrt( 1 - cosTh * cosTh );
        
        if (cosTh > 0.) {
            zInc = patchTop;
        }
        else {
            zInc = patchBottom;
        }
        xInc = patchXmin + patchWidX * RandFlat::shoot();
        yInc = patchYmin + patchWidY * RandFlat::shoot();
    }
    
    Point IncPoint(xInc, yInc, zInc);
    m_launchPoint = IncPoint;
    m_launchDir = Vector( -(cos(phi) * sinTh), -(sin(phi) * sinTh), -(cosTh));
}

void FluxSource::getGalacticDir(double l,double b){
    
    //here is the new mechanism:
    double theta=sqrt(pow(b,2)*pow(l,2));
    if (theta==0.){theta+=0.000000000001;}  //to fix divide-by-zero errors
    double phi=acos(l/theta);
    //here we construct the cartesian galactic matrix
    Vector gamgal(sin(theta)*cos(phi) , sin(theta)*sin(phi) , cos(theta));

    //get the transformation matrix..
    Rotation galtoglast=GPS::instance()->orbit()->CELtransform(GPS::instance()->time());

    //and do the transform:
    setLaunch(galtoglast*gamgal);
  
    /*  old way to do it - the galToGlast function is depreciated
    std::pair<double,double> v;
    v=GPS::instance()->galToGlast(std::make_pair<double,double>(l,b));
    
    double theta=v.first;
    double phi=v.second;
    //std::cout << "getting galactic direction..." << std::endl;
    setLaunch(theta,phi);
    */
    m_launch = GALACTIC;
}


std::string FluxSource::title () const
{
    std::strstream t;
    t << m_spectrum->title() << ", ";
    switch (m_launch) {
    case NONE:
        t << "range(" << _minCos << ',' << _maxCos << "), ";
        if( theta() == 0) break;
    case DIRECTION:
        t << "angle(" << theta()*180/M_PI << ','
            << phi()*180/M_PI << "), ";
        break;
    case POINT:
        t << "launch(" << m_launchDir << m_launchPoint << "), ";
        break;
    case SURFACE:
        t << "box(" << patchXmin << ',' << patchXmax << ',' << patchYmin << ',' 
            << patchYmax << ',' << patchTop << ',' << patchBottom << ") theta(" 
            << _minCos << ',' << _maxCos << "), phi(" << _minPhi << ',' 
            << _maxPhi << "), ";
        break;
    case SPECTRUM:
        break;
    case GALACTIC:
        break;
    case PATCHFIXED:
        t << "box(" << patchXmin << ',' << patchXmax << ',' << patchYmin 
            << ',' << patchYmax << ',' << patchTop << ',' << patchBottom 
            << ") theta(" << _theta << "), phi(" << _phi << "), ";
    }
    t << "flux("<< flux(0.) << ")" << '\0';    
    std::string s(t.str());
#ifdef WIN32
    t.freeze(false);
#endif
    return s;
}

/*
void FluxSource::transformDirection(){

 switch (m_frametype) {
    case GLAST:
        m_transformDir=m_launchDir;
        break;
    case EARTH:
        m_transformDir=GPS::instance()->earthToGlast(m_launchDir);
        break;
    case GALAXY:
        m_transformDir=GPS::instance()->galaxyToGlast(m_launchDir);
        break;
 }
 
}*/

void FluxSource::refLaunch(LaunchType launch) {m_launch=launch;}
void FluxSource::refPoint(PointType point) {m_pointtype=point;}

double FluxSource::interval (double time){
return m_spectrum->interval(time);
}

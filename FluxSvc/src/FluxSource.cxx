//$Header$

#include "FluxSvc/FluxSource.h"

#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "xml/Dom.h"
#include "CLHEP/Random/RandFlat.h"

#include "astro/SkyDir.h"

#include "SpectrumFactoryTable.h"
#include "SimpleSpectrum.h"

#include "FluxException.h" // for FATAL_MACRO
#include "GPS.h"

#include <algorithm>
#include <sstream>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class LaunchPoint
    @brief nested launch strategy base class for point determination
    
    The virtual base class manages the point itself
*/
class FluxSource::LaunchPoint  { 
public:
    LaunchPoint(){}
    LaunchPoint(const HepPoint3D& pt):m_pt(pt){}
    virtual ~LaunchPoint(){}

    /// access to direction, perhaps set by the execute()
    virtual const HepPoint3D& point()const {return m_pt;}
    const HepPoint3D& operator()()const{return point();}

    /// execute the strategy, perhaps depending on direction
    virtual void execute(const HepVector3D& dir){};

    /// set the point
    void setPoint(const HepPoint3D& pt){ m_pt = pt;}

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << "point" << m_pt;
        return t.str();
    }

private:
    HepPoint3D m_pt;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomPoint
    @brief nested launch strategy derived class
    This is the standard strategy, which takes a direction and creates a point in
    a disk centered at the origin with area 6 m^2 (or so)
*/
class FluxSource::RandomPoint : public LaunchPoint{ 
public:
    RandomPoint(double radius, double backoff)
        :m_radius(radius), m_backoff(backoff)
    { 

    }

    virtual void execute(const HepVector3D& dir){
        HepRotation r_pln;

        double ly = dir.y(), lx = dir.x();
        if( lx !=0 || ly !=0 ) { 
            r_pln.rotate(acos(dir.z()),  HepVector3D(-ly, lx, 0.));
        }

        // pick a random position on the planar section of a sphere through 
        // its midpoint
        double 
            azimuth = RandFlat::shoot( 2*M_PI ),
            rad = m_radius*(sqrt(RandFlat::shoot()));

        // create two vectors to describe the particle launch: one to describe
        // the point in the plane perpendicular to the launch direction (within
        // the cross-section of the sphere containing the instrument) and a 
        // second to describe the distance along the normal between the launch 
        // point and that plane.
        HepPoint3D posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

        // define actual launch point
        setPoint( r_pln*posLaunch - m_backoff*dir);
    }

    /// return info, 
    std::string title()const{
        std::stringstream t;
        t << "radius("<< m_radius << ")";
        return t.str();
    }


private:
    double m_radius;
    double m_backoff;

}; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class FixedPoint
    @brief nested launch strategy derived class
    
    This strategy uses a fixed launch point for a pencil beam. If the radius is nonzero,
    the beam will be spread out uniformly on a disk perpendicular to the incoming direction
*/
class FluxSource::FixedPoint : public LaunchPoint{ 
public:
    FixedPoint( const HepPoint3D& pt, double radius)
        :  LaunchPoint(pt)
        ,  m_disk_radius(radius)
        ,  m_base_point(pt)
    {}

    virtual void execute(const HepVector3D& dir){
        if(m_disk_radius==0) return; // just use 
    
        HepRotation r_pln;

        double ly = dir.y(), lx = dir.x();
        if( lx !=0 || ly !=0 ) { 
            r_pln.rotate(acos(dir.z()), HepVector3D(-ly, lx, 0.));
        }
        double 
            azimuth = RandFlat::shoot( 2*M_PI ),
            rad = m_disk_radius*(sqrt(RandFlat::shoot()));
        HepPoint3D posLaunch(rad*cos(azimuth), rad*sin(azimuth), 0.);

        setPoint(r_pln*posLaunch + m_base_point);

    };
    virtual std::string title() const {
        if( m_disk_radius==0) return LaunchPoint::title();
        std::stringstream t;
        t << ", radius(" << m_disk_radius << ")";
        return LaunchPoint::title() + t.str();
    }
private:
    double m_disk_radius;
    HepPoint3D m_base_point;
};  
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Patch
    @brief nested launch strategy derived class
    Gets a point randomly from a box
*/
class FluxSource::Patch : public FluxSource::LaunchPoint{ 
public:
    Patch( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
       :m_xmin(xmin), m_dx(xmax-xmin), 
        m_ymin(ymin), m_dy(ymax-ymin), 
        m_zmin(zmin), m_dz(zmax-zmin)
    {
    }

    virtual void execute(const HepVector3D& ){
        setPoint(HepPoint3D( 
            m_xmin + m_dx*RandFlat::shoot(),
            m_ymin + m_dy*RandFlat::shoot(),
            m_zmin + m_dz*RandFlat::shoot()) );
    }
    virtual std::string title() const {
        std::stringstream t;
        t << "patch(" 
            << m_xmin << "," << m_xmin+m_dx << ","
            << m_ymin << "," << m_ymin+m_dy << ","
            << m_zmin << "," << m_zmin+m_dz << ")" ;
        return t.str();
    }

private:
    double m_xmin, m_dx, m_ymin, m_dy, m_zmin, m_dz;    
}; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class LaunchDirection
    @brief nested launch strategy base class
*/
class FluxSource::LaunchDirection  {
public:
    LaunchDirection():m_skydir(false){}


    LaunchDirection(double theta, double phi)
        :m_skydir(false)
    {
        HepVector3D dir(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
        setDir(-dir); // minus due to z axis pointing UP!
    }
    LaunchDirection(astro::SkyDir sky)
        :m_skydir(true)
    {
        m_dir = sky.dir();
    }
    /** @brief choose a direction
        @param KE kinetic energy
        @parem time mission time
        */
    virtual void execute(double KE, double time){
        //TODO: account for transformation?
    }

    const HepVector3D& operator()()const {return dir();}

    virtual const HepVector3D& dir()const { return m_dir;}
    void setDir(const HepVector3D& dir){m_dir=dir;}


    //! solid angle: default of 1. for a point source
    virtual double solidAngle()const {
        return 1.;
    }

    /// return info, default if not overriden
    virtual std::string title()const{
        std::stringstream t;
        t << " dir" << m_dir ;
        return t.str();
    }

private:
    HepVector3D m_dir;
    bool  m_skydir;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class RandomDirection
    @brief nested launch strategy derived class
    Gets a point randomly from a box
*/
class FluxSource::RandomDirection : public FluxSource::LaunchDirection{ 
public:
    /** ctor:
    @param minc  minimum value of cos(theta)
    @param maxc  maximum value of cos(theta)
    @param theta 
    */
    RandomDirection(double minc, double maxc, double theta=0, double phi=0)
        : m_theta(theta)
        , m_phi(phi)
    {
    using std::min;
    using std::max;

        // require _maxCos > _minCos for solid angle calculation
        m_minCos   = min(minc,maxc);
        m_maxCos   = max(minc,maxc);
        if(m_minCos==m_maxCos) {
            if(m_minCos!=-1) m_minCos-=0.001; else m_maxCos +=0.001;
        }
        m_minPhi = 0; m_maxPhi=2*M_PI;

    }

    virtual void execute(double, double){
            double  costh = -RandFlat::shoot(m_minCos, m_maxCos),
                    sinth = sqrt(1.-costh*costh),
                    phi = RandFlat::shoot(m_minPhi, m_maxPhi);
            
            HepVector3D dir(cos(phi)*sinth, sin(phi)*sinth, costh);

           // extra rotation in case not zenith pointing (beware, might be
            // confusing)
            // keep x-axis perpendicular to zenith direction
            if (m_theta != 0.0) dir.rotateX(m_theta).rotateZ(m_phi);
            setDir(dir);

    }
       //! solid angle
    virtual double solidAngle()const {
        return 2*M_PI*(m_maxCos-m_minCos);
    }


    virtual std::string title()const {
        std::stringstream t;
        t << "range(" << m_minCos << ',' << m_maxCos << ") ";
        if( m_theta != 0){
            t << ", angle(" << m_theta*180/M_PI << ',' << m_phi*180/M_PI << ") ";
        }
        return t.str();
    }


private:
    double m_minCos, m_maxCos;
    double m_minPhi, m_maxPhi;
    double m_theta, m_phi;
    
}; 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class SourceDirection
    @brief nested launch strategy derived class
    Gets a direction from the ISpectrum class
*/
class FluxSource::SourceDirection : public FluxSource::LaunchDirection{ 
public:
    /** Ctor:
    @param spectrum pointer to the ISpectrum object that will provide the direction
    @param galactic if true, interpret pair as l,b (in degrees); otherwise costh, phi
    */
    SourceDirection(ISpectrum* spectrum, bool galactic)
        : m_spectrum(spectrum)
        , m_galactic(galactic){}

        void execute(double ke, double time){


            std::pair<float,float> direction 
                    = m_spectrum->dir(ke,HepRandom::getTheEngine());

            if( !m_galactic) {
                // special option that gets direction from the spectrum object
                // note extra - sign since direction corresponds to *from*, not *to*

                double  costh = direction.first,
                        sinth = sqrt(1.-costh*costh),
                        phi = direction.second;
                setDir(-HepVector3D(cos(phi)*sinth, sin(phi)*sinth, costh));

            }else {
                // iterpret direction as l,b for a galactic source
                double  l = direction.first,
                        b = direction.second;
                //then set up this direction:
                astro::SkyDir unrotated(l,b,astro::SkyDir::GALACTIC);
                //get the transformation matrix..
                HepRotation celtoglast
                    =GPS::instance()->transformCelToGlast(GPS::instance()->time() + time);

                //and do the transform:
                setDir(celtoglast*unrotated());
            }
        }

        //! solid angle
        virtual double solidAngle()const {
            return m_spectrum->solidAngle();
        }
        
        virtual std::string title()const {
            return "(use_spectrum)";
        }

private:
    ISpectrum* m_spectrum;
    bool   m_galactic;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                   FluxSource constructor
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FluxSource::FluxSource(const DOM_Element& xelem )
: EventSource (xelem)
, m_spectrum(0)
{
    static double d2r = M_PI/180.;


    ISpectrum*   s = 0;
    std::string class_name;
    std::string source_params; 
    
    DOM_Element   spec = xml::Dom::findFirstChildByName(xelem, "spectrum");
    
    if (spec == DOM_Element()) {
        
        // source has no imbedded spectrum element: expect a name
        class_name = xml::Dom::transToChar(xelem.getAttribute("name"));
#if 0
        useSpectrumDirection(); // and will use its direction generation
#else
//???        m_launch_dir = new SourceDirection(m_spectrum, frame=="galaxy"); 

#endif   
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
        std::string spectrum_name = xml::Dom::transToChar(spec.getAttribute("name"));
        std::string spectrum_energyscale = xml::Dom::transToChar(spec.getAttribute("escale"));
        
        if(spectrum_energyscale == "GeV"){ m_energyscale=GeV;
        }else if(spectrum_energyscale == "MeV"){ m_energyscale=MeV;
        }else{
            std::cout << "bad energy scale declaration on spectrum:"
                << spectrum_energyscale << " , exiting.";
            return;} //this line "just in case"
        
        
        if (typeTagName.equals("particle")) s = new SimpleSpectrum(specType);
        else if (typeTagName.equals("SpectrumClass")) {
            // attribute "name" is the class name
            class_name = xml::Dom::transToChar(specType.getAttribute("name"));
            source_params= xml::Dom::transToChar(specType.getAttribute("params"));
        }
        else {
            // no, the tag itself
            class_name = xml::Dom::transToChar(typeTagName);//.transcode();
        }

        //if s is still 0, we need to create the internal spectrum object.
        if( s==0) {
            //		std::vector<float> paramvec; parseParamList(source_params, paramvec);
            s = SpectrumFactoryTable::instance()->instantiate(class_name, source_params);
            if(s==0){
                
                std::cerr << "List of known Spectrum classes:\n" ;
                std::list<std::string>list= SpectrumFactoryTable::instance()->spectrumList();
                for( std::list<std::string>::iterator i = list.begin(); i!=list.end(); ++i)
                    std::cerr << "\t" << *i << std::endl;
                FATAL_MACRO("Unknown Spectrum: "<< class_name);
                return;
            }
        }
        m_spectrum =s;


        // second child element is angle
        DOM_Element angles = xml::Dom::getSiblingElement(specType);
        DOMString anglesTag = angles.getTagName();
        
        if (anglesTag.equals("solid_angle") ) {
            m_launch_dir = new RandomDirection(
                atof(xml::Dom::transToChar(angles.getAttribute("mincos"))),
                atof(xml::Dom::transToChar(angles.getAttribute("maxcos"))),
                atof(xml::Dom::transToChar(angles.getAttribute("theta"))) * d2r, 
                atof(xml::Dom::transToChar(angles.getAttribute("phi"))) *d2r);

        }
        else if (anglesTag.equals("direction") ) {
            m_launch_dir = new LaunchDirection(
                atof(xml::Dom::transToChar(angles.getAttribute("theta"))) * d2r, 
                atof(xml::Dom::transToChar(angles.getAttribute("phi"))) *d2r);
        }
        else if (anglesTag.equals("use_spectrum") ) {
            std::string frame = xml::Dom::transToChar(angles.getAttribute("frame"));
            m_launch_dir = new SourceDirection(m_spectrum, frame=="galaxy"); 
        }
        else if(anglesTag.equals("galactic_dir")){
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                    atof(xml::Dom::transToChar(angles.getAttribute("l"))),
                    atof(xml::Dom::transToChar(angles.getAttribute("b"))), 
                    astro::SkyDir::GALACTIC
                    ) 
                );
        }
        else if(anglesTag.equals("celestial_dir")){
            m_launch_dir = new LaunchDirection(
                astro::SkyDir(
                    atof(xml::Dom::transToChar(angles.getAttribute("ra"))),
                    atof(xml::Dom::transToChar(angles.getAttribute("dec"))), 
                    astro::SkyDir::CELESTIAL 
                    ) 
                );

        }
        else if(anglesTag.equals("galactic_spread")){
            FATAL_MACRO("not implemented");
        }
        else {
            FATAL_MACRO("Unknown angle specification in Flux::Flux \""
                << xml::Dom::transToChar(anglesTag) << "\"" );
        }
        
        // third child element is optional launch spec
        DOM_Element launch = xml::Dom::getSiblingElement(angles);
        
        if(launch !=DOM_Element()) {
            DOMString launchTag = launch.getTagName();
            
             if(launchTag.equals("launch_point")){
                m_launch_pt = new FixedPoint(
                    HepPoint3D(
                        atof(xml::Dom::transToChar(launch.getAttribute("x"))),
                        atof(xml::Dom::transToChar(launch.getAttribute("y"))),
                        atof(xml::Dom::transToChar(launch.getAttribute("z"))) ),
                    atof(xml::Dom::transToChar(launch.getAttribute("beam_radius")))
                    );

                
            }else if(launchTag.equals("patch")){
                m_launch_pt = new Patch( 
                    atof(xml::Dom::transToChar(launch.getAttribute("xmax"))),
                    atof(xml::Dom::transToChar(launch.getAttribute("xmin"))),
                    atof(xml::Dom::transToChar(launch.getAttribute("ymax"))),
                    atof(xml::Dom::transToChar(launch.getAttribute("ymin"))), 
                    atof(xml::Dom::transToChar(launch.getAttribute("zmax"))),
                    atof(xml::Dom::transToChar(launch.getAttribute("zmin"))) 
                    );
                
            }else {
                FATAL_MACRO("Unknown launch specification in Flux::Flux \""
                    << xml::Dom::transToChar(launchTag) << "\"" );
            }
        } else {
            double radius = sqrt(totalArea() / M_PI ) * 1000;    // radius in mm
            m_launch_pt = new RandomPoint(radius, radius);
        }
    }
}


FluxSource::~FluxSource()
{
    delete m_spectrum;
    delete m_launch_pt;
    delete m_launch_dir;

}

void FluxSource::spectrum(ISpectrum* s, double emax)
{
    if (emax > 0) {
        //m_rmax =  s->fraction(emax);
        std::cerr << "exercising obsolete function fraction" << std::endl;
    }
    m_spectrum = s;
    //const char* name = s->particleName();
}

FluxSource* FluxSource::event(double time)
{
    // Purpose and Method: generate a new incoming particle
    // Inputs  - current time
    // Outputs - pointer to the "current" fluxSource object.
    m_interval = calculateInterval(time);
    computeLaunch(time+m_interval);

    //now set the actual interval to be what FluxMgr will get
    EventSource::setTime(time+m_interval);

    // Finally rock it.
    HepRotation correctForTilt =GPS::instance()->rockingAngleTransform(GPS::instance()->time());
    m_correctedDir = correctForTilt*m_launchDir;

    return this;
}

double FluxSource::calculateInterval (double time)
{   
    double interval=m_spectrum->interval(time);
    if( interval>0 ){
        // the spectum computed an interval: use it
        return interval;
    }

    // otherwise do a Poison from the the flux, solid angle, and area factor
    return -log(1.-RandFlat::shoot(1.))/rate(time);
}

double FluxSource::flux(double time) const
{
    if(!enabled()){ return 0;}
    if(m_spectrum->flux(time)){ return m_spectrum->flux(time);}
    else{return EventSource::flux(time);}
}

double FluxSource::rate(double time) const
{
    //TODO: area is only relevant for RandomPoint, and flux per unit area
    return m_launch_dir->solidAngle() * flux(time) * totalArea();
}

void FluxSource::computeLaunch (double time)
{
    // get the KE from the spectrum object
    m_energy = spectrum()->energySrc( HepRandom::getTheEngine(), time );

    // convert to MeV if necessary
    if(m_energyscale==GeV){
        m_energy *= 1000.;
    }
    
    // set launch direction , position (perhaps depending on direction)
    m_launch_dir->execute(m_energy, time);
    m_launch_pt->execute((*m_launch_dir)());
    
    m_launchPoint = (*m_launch_pt)();
    m_launchDir  = (*m_launch_dir)();

    //correctForTiltAngle();
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


std::string FluxSource::title () const
{
    return m_spectrum->title() + ", "
        +  m_launch_pt->title() +", "
        +  m_launch_dir->title();
}


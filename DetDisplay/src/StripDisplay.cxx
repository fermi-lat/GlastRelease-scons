/** 
* @file StripDisplay.cxx
* @brief Declaration and definition of the algorithm StripDisplay.
*
*  $Header$
*/
// Original Author: T. Burnett

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/GaudiException.h" 

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/MonteCarlo/McTkrStrip.h"

//gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"
#include "geometry/CoordTransform.h"
#include "geomrep/BoxRep.h"
#include "geometry/Point.h"
#include "GuiSvc/IGuiTool.h"

#include <cassert>

/*! \class StripDisplay
\brief Display of strip hits.

*/
namespace {
    // local pointers to GlastDetSvc and TkrGeometrySvc
    IGlastDetSvc* _detSvc;
    ITkrGeometrySvc* _tkrGeoSvc;
}

class StripDisplay :
    public AlgTool, 
    virtual public IIncidentListener, 
    virtual public IGuiTool 
{

public:
    StripDisplay ( const std::string& type, const std::string& name, const IInterface* parent);
    ~StripDisplay();

    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    

    // these needed to implement gui::MenuClient interface 
    void quit(){};
    void finishSetup(){};  // dummy

    class StripRep;
private:
    void clear();
    void displayStrips(const Event::McTkrStripCol& strips);

    /// manage incidents
    void beginEvent();
    void endEvent();
    IDataProviderSvc* eventSvc() const {
        if ( 0 == m_EDS ) {
            StatusCode sc = service( "EventDataSvc", m_EDS, true );
            if( sc.isFailure() ) {
                throw GaudiException("Service [EventDataSvc] not found", name(), sc);
            }
        }
        return m_EDS;
    }     
    mutable IDataProviderSvc* m_EDS;              ///< Event data service

    StripRep* m_stripRep;
    StripRep* m_noiseRep;
    /// Strip data declaration
    class Strip 
    {
    public:
        // constructor
        Strip (int index = -1, double energy = 0, bool noise=false,
            double deltaX=0., double deltaY=0.)
            :m_index(index), m_energy(energy), m_noise(noise),
            m_deltaX(deltaX), m_deltaY(deltaY)
        {}

        // access, 
        void addEnergy (float e)    {  m_energy += e;   }   // add energy
        void energy (float e)       {  m_energy = e;    }   // set energy
        float energy () const       {  return m_energy; }   // get energy
        unsigned int index() const  {  return m_index;  }   // get index
        bool noise() const          {  return m_noise;  }   // get noise status
        double getDeltaX() const    {  return m_deltaX; }
        double getDeltaY() const    {  return m_deltaY; }

        // static parameters
        // undefined strip (non-existent)
        static unsigned int undef_strip () { return 65535; } 

    private:
        int     m_index;  // strip number, -1 if invalid
        float   m_energy; // charge deposited
        bool    m_noise;  // true if generated by noise
        double  m_deltaX;  
        double  m_deltaY;
    };

    typedef  std::vector<Strip>  SiDetector;
    typedef std::map<idents::VolumeIdentifier, SiDetector > SiDetMap;
    SiDetMap m_SiMap;
};

static const ToolFactory<StripDisplay>  Factory;
const IToolFactory& StripDisplayFactory = Factory;

StripDisplay::StripDisplay(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent)
                           : AlgTool( type, name, parent ) 
                           ,m_EDS(0) 
{
    // declare the properties

    // Declare additional interface
    declareInterface<IGuiTool>(this);
}

//   Strip display rep
class   StripDisplay::StripRep : public gui::DisplayRep {
public:
    StripRep(bool noise):m_count(0), m_noise(noise){}
    void   accept ( const SiDetector& d, const CoordTransform& T_plane )
    {
        // static variable implmentation
        static float pScale = 50.f;

        if ( d.empty()) return;
        // Now draw the hits as red thermometers bars
        setColor(m_noise?  "orange": "red");   

        for( SiDetector::const_iterator strip = d.begin(); strip != d.end(); ++strip) {
            if (strip->noise() == m_noise) // select noise or not
            {
                float h = strip->energy() * pScale;
                float x0 = _detSvc->stripLocalX(strip->index()) + strip->getDeltaX();
                float y0 = _tkrGeoSvc->trayWidth()/2 + strip->getDeltaY();
                Point from(x0,y0,0), to(x0,y0,h);
                from.transform(T_plane); to.transform(T_plane);
                moveTo(from);
                lineTo(to);
            }
        }

        setColor("black");
        ++m_count;
    }
    void clear(){
        m_count=0;
        DisplayRep::clear();
    }
    void update(){};
    int count()const { return m_count;}
private:
    int m_count;
    bool m_noise;
};

StripDisplay::~StripDisplay()
{
}

StatusCode StripDisplay::initialize(gui::GuiMgr* guiMgr) 
{
    using namespace gui;   
    MsgStream log(msgSvc(), name());

    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    StatusCode status = service ("IncidentSvc", incsvc, true);
    if( status.isFailure() ) {
        throw GaudiException("Service [IncidentSvc] not found", name(), status);
    }
    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    // now try to find the GlastDevSvc service
    StatusCode sc = service("GlastDetSvc", _detSvc, true);

    if (!sc.isSuccess ()){
        throw GaudiException("Service [GlastDetSvc] not found", name(), status);
    }
    // now try to find the TkrGeometrySvc service
    sc = service("TkrGeometrySvc", _tkrGeoSvc, true);

    if (!sc.isSuccess ()){
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), status);
    }

    gui::DisplayControl& dc = guiMgr->display();

    gui::DisplayControl::DisplaySubMenu& SiHitsMenu = dc.subMenu("Si Strip Hits");
    SiHitsMenu.add(m_stripRep = new StripRep(false), "Strip hits");
    SiHitsMenu.add(m_noiseRep = new StripRep(true), "Noise hits");
    return sc;
}

void StripDisplay::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
}

void StripDisplay::beginEvent()
{
    // clear the list of hit detectors
    clear();
}

void StripDisplay::endEvent()
{
    using namespace Event;

    SmartDataPtr<Event::McTkrStripCol> stripCol(eventSvc(), EventModel::MC::McTkrStripCol);
    if( stripCol ){
        displayStrips (stripCol);
    }
}

void StripDisplay::clear()
{ 
    for( SiDetMap::iterator im1 = m_SiMap.begin(); im1 != m_SiMap.end(); ++im1){
        im1->second.clear();
    }
}

void StripDisplay::displayStrips(const Event::McTkrStripCol& strips)
{
    MsgStream   log( msgSvc(), name() );
    using namespace Event;

    // first sort into planes according to ID
    for(McTkrStripCol::const_iterator ihit=strips.begin(); ihit != strips.end(); ++ihit) {
        const McTkrStrip & strip = **ihit;
        m_SiMap[strip.getId()].push_back(Strip(strip.getStripNumber(),strip.energy(),
            strip.getNoise(), strip.getDeltaX(), strip.getDeltaY()));
    }

    // now loop over planes; get transformation for each one,  pass info to the display rep.
    for(SiDetMap::iterator im = m_SiMap.begin(); im != m_SiMap.end(); ++im){
        const SiDetector& strips = im->second;
        if (strips.empty()) continue;
        idents::VolumeIdentifier plane_id= im->first;

        //  form a plane tranformation by averaging the offsers for  
        //       wafers (0,0) and (nWafers-1,nWafers-1)
        //  bit of a kludge, but it works!
        idents::VolumeIdentifier   waferMin_id=plane_id, waferMax_id=plane_id;
        waferMin_id.append(0); waferMin_id.append(0);
        int topWafer = _tkrGeoSvc->nWaferAcross() - 1;
        waferMax_id.append(topWafer); waferMax_id.append(topWafer);
        HepTransform3D TMin, TMax;
        _detSvc->getTransform3DByID(waferMin_id, &TMin);
        _detSvc->getTransform3DByID(waferMax_id, &TMax);
        Hep3Vector offset = 0.5*(TMin.getTranslation() + TMax.getTranslation());
        CoordTransform  T_plane(TMin.getRotation(), offset);

        // pass the strip list and transformation to the two instances of the rep
        m_stripRep->accept(strips, T_plane);
        m_noiseRep->accept(strips, T_plane);
    }
}
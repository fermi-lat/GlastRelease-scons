/** @file MCdisplay.cxx 
     @brief Declatation, definition of class MCdisplay
  $Header$
*/

// includes
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McRelTableDefs.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "idents/VolumeIdentifier.h"

#include "GuiSvc/IGuiTool.h"

#include "gui/DisplayControl.h"
#include "gui/SubMenu.h"
#include "gui/GuiMgr.h"
#include "gui/Command.h"

#include "geometry/Box.h"
#include "geometry/CoordTransform.h"
#include "geomrep/BoxRep.h"

#include "CLHEP/Geometry/Transform3D.h"
#include <map>


/** 
* @class MCdisplay
*
* @brief  A Tool that create a display of  Monte Carlo items on the TDS 
* $Header$
*/
class MCdisplay : public AlgTool, 
    virtual public IIncidentListener, 
    virtual public IGuiTool 
{
public:

    //! Create the tool
    //! @param name The name of the service
    //!
    MCdisplay ( const std::string& type, const std::string& name, const IInterface* parent);

    virtual ~MCdisplay (){};

    /// implement to define gui elements: will be called from GuiSvc
    StatusCode initialize(gui::GuiMgr*);

    /// Handles incidents, implementing IIncidentListener interface
    virtual void handle(const Incident& inc);    

    // these needed to implement gui::MenuClient interface 
    void quit(){};
    void finishSetup(){};  // dummy
private:

    //! add a hit detector box to the display
    void displayBox(idents::VolumeIdentifier id, std::string category);
    //! add a hit to the display 
    //! @param a,b initial, final points for the step
    void addHit(const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b);

    //! Add rep to display the tracks
    //! @param track a list of dots to connect
    //! @param charge tell if neurtal
    //! @param primary flag if primary track, to distinguish from neutrals?
    typedef std::vector<Event::McTrajectoryPoint*> PointList;
    void addTrack(const PointList& track,  int charge, bool primary=false);

    gui::DisplayControl* m_display;

    //! access all DisplayRep pointers by name in this map
    std::map<std::string, gui::DisplayRep*> m_detmap;

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
    mutable IGlastDetSvc* m_gdsvc;
    mutable IParticlePropertySvc* m_ppsvc;

    typedef std::map<idents::VolumeIdentifier, unsigned int> DetectorList;

  
    /// A map to keep track of hit detectors for display
    DetectorList m_detectorList;
#if defined(__GNUC__) && (__GNUC__ == 2)
    //! minimal rep just to append stuff to. Allows default color
    class EmptyRep : public gui::DisplayRep {
    public:
        EmptyRep(std::string color="black"):m_color(color){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor(m_color);}
    private: 
        std::string m_color;
    };

    class TrackRep : public gui::DisplayRep {
    public:
        TrackRep( const MCdisplay::PointList& track, bool neutral=false, bool primary=false){
            if( neutral) setColor("white");
            MCdisplay::PointList::const_iterator pit = track.begin();
            if(primary) markerAt(*pit);
            moveTo(*pit++);
            for(; pit !=track.end(); ++pit) lineTo(*pit);
        }
        void update(){}
    };

    class LineRep : public gui::DisplayRep {
    public:
        LineRep(const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b) 
        {
            markerAt(a); moveTo(a);  lineTo(b);
        }
        void update(){}
    };
#endif

};


// declare the service factories for the MCdisplay
static ToolFactory<MCdisplay> a_factory;
const IToolFactory& MCdisplayFactory = a_factory;

//____________________________________________________________________________
///       Standard Constructor
MCdisplay::MCdisplay(const std::string& type, 
                     const std::string& name, 
                     const IInterface* parent)
                     : AlgTool( type, name, parent ) 
                     ,m_EDS(0) , m_gdsvc(0), m_ppsvc(0)
{
    // declare the properties

    // Declare additional interface
    declareInterface<IGuiTool>(this);
}


//____________________________________________________________________________
// initialize
StatusCode MCdisplay::initialize (gui::GuiMgr* guiMgr) 
{
    StatusCode  status = StatusCode::SUCCESS;

    // bind all of the properties for this service
    status = setProperties ();

    // open the message log
    MsgStream log( msgSvc(), name() );
    gui::DisplayControl& dc = guiMgr->display();
    gui::DisplayControl::DisplaySubMenu& m = dc.subMenu("MonteCarlo");
#if !defined(__GNUC__) || (__GNUC__ != 2)
    //! minimal rep just to append stuff to. Allows default color
    class EmptyRep : public gui::DisplayRep {
    public:
        EmptyRep(std::string color="black"):m_color(color){}
        void update(){}
        void clear(){DisplayRep::clear(); setColor(m_color);}
    private: 
        std::string m_color;
    };
#endif
    m.add(m_detmap["steps"] = new EmptyRep("red"), "hits", false);

    m.add(m_detmap["hit_boxes"] = new EmptyRep("orange"), 
        "hit Pos detectors");

    m.add(m_detmap["integrating_boxes"] = new EmptyRep("blue"), 
        "hit Int detectors");

    m.add((m_detmap["primary"] = new EmptyRep("black")), "primary track");
    m.add((m_detmap["tracks"]= new EmptyRep("black")), "charged tracks");
    m.add((m_detmap["neutrals"] = new EmptyRep("white")), "neutral tracks");

    // use the incident service to register begin, end events
    IIncidentSvc* incsvc = 0;
    status = service ("IncidentSvc", incsvc, true);
    if( status.isFailure() ) {
        throw GaudiException("Service [IncidentSvc] not found", name(), status);
    }
    incsvc->addListener(this, "BeginEvent", 100);
    incsvc->addListener(this, "EndEvent", 0);

    if( (status=service("ParticlePropertySvc", m_ppsvc)).isFailure() ) {
        throw GaudiException("Service [ParticlePropertySvc] not found", name(), status);
    }

    // now try to find the GlastDevSvc service
    status = service("GlastDetSvc", m_gdsvc);
    if( status.isFailure() ) {
        throw GaudiException("Service [GlastDetSvc] not found", name(), status);
    }
    return status;
}
// handle "incidents"
//____________________________________________________________________________
void MCdisplay::handle(const Incident &inc)
{
    if( inc.type()=="BeginEvent")beginEvent();
    else if(inc.type()=="EndEvent")endEvent();
}
//____________________________________________________________________________
void MCdisplay::displayBox(idents::VolumeIdentifier id, std::string category)
{
    if(m_detectorList[id]>0)return; // already displayed
    ++m_detectorList[id];
    HepGeom::Transform3D T;
    m_gdsvc->getTransform3DByID(id, &T);

    std::string shape;
    std::vector<double> params;
    m_gdsvc->getShapeByID(id, &shape, &params); 
    if(shape != "box") return;
    Box b( params[0], params[1], params[2]);
    b.transform(CoordTransform(T.getRotation(), T.getTranslation()));
    m_detmap[category]->append(BoxRep(b));

}
//____________________________________________________________________________
void MCdisplay::addTrack(const PointList & track, int charge, bool primary)
{
#if !defined(__GNUC__) || (__GNUC__ !=2)
    class TrackRep : public gui::DisplayRep {
    public:
        TrackRep( const MCdisplay::PointList& track, bool neutral=false, bool primary=false){
            if( neutral) setColor("white");
            MCdisplay::PointList::const_iterator pit = track.begin();
            CLHEP::Hep3Vector hit = (*pit)->getPoint();
            if(primary) markerAt(hit);
            moveTo(hit);
            pit++;
            for(; pit !=track.end(); ++pit) lineTo((*pit)->getPoint());
        }
        void update(){}
    };
#endif
    if( primary)        m_detmap["primary"]->append(TrackRep(track, charge==0, primary));
    else  if(charge==0) m_detmap["neutrals"]->append(TrackRep(track));
    else                m_detmap["tracks"]->append(TrackRep(track));
}
//____________________________________________________________________________
void MCdisplay::addHit( const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b)
{
#if !defined(__GNUC__) || (__GNUC__ !=2)
    class LineRep : public gui::DisplayRep {
    public:
        LineRep(const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b) 
        {
            markerAt(a); moveTo(a);  lineTo(b);
        }
        void update(){}
    };
#endif
    m_detmap["steps"]->append( LineRep(a,b) );
}

//____________________________________________________________________________
void MCdisplay::beginEvent()
{
  // clear the list of hit detectors
  m_detectorList.clear();
}
//____________________________________________________________________________
void MCdisplay::endEvent()
{
    SmartDataPtr<Event::EventHeader>   header(eventSvc(),    EventModel::EventHeader);
    if(! header) return; // should not happen
    SmartDataPtr<Event::MCEvent>     mcheader(eventSvc(),    EventModel::MC::Event);

    int mc_src_id = mcheader->getSourceId();
    int mcevent = mcheader->getSequence();
    double time = header->time();
    int event = header->event();
    MsgStream log( msgSvc(), name() );

    // ---------- position hits -----------------------------------------
    SmartDataPtr<Event::McPositionHitVector>  posHits(eventSvc(), "/Event/MC/PositionHitsCol");
    if( posHits)   for(Event::McPositionHitVector::const_iterator ihit=posHits->begin();
        ihit != posHits->end(); ihit++)
    {
        idents::VolumeIdentifier id = (*ihit)->volumeID();
        displayBox(id, "hit_boxes");
        // add the hit
        HepGeom::Transform3D T;
        m_gdsvc->getTransform3DByID(id, &T);
        HepPoint3D entry = T.getTranslation() + T.getRotation()*(*ihit)->entryPoint();
        HepPoint3D exit = T.getTranslation() + T.getRotation()*(*ihit)->exitPoint();
        addHit(entry, exit);

    }      
    // ---------- integrating  hits -----------------------------------------
    SmartDataPtr<Event::McIntegratingHitVector> intHits(eventSvc(), "/Event/MC/IntegratingHitsCol");

   if(intHits) for(Event::McIntegratingHitVector::const_iterator inHit=intHits->begin(); 
        inHit != intHits->end(); inHit++) 
    {
        idents::VolumeIdentifier id = (*inHit)->volumeID();
        displayBox(id, "integrating_boxes" );
    }
    // ---------- tracks -----------------------------------------
    bool primary=true;

    // If there are trajectories in the TDS, we use them       
    SmartDataPtr<Event::McPartToTrajectoryTabList> 
        mcPartToTraj(eventSvc(), "/Event/MC/McPartToTrajectory");

    if (mcPartToTraj != 0)
    {
        for(Event::McPartToTrajectoryTabList::const_iterator trajIter = mcPartToTraj->begin(); 
            trajIter != mcPartToTraj->end(); trajIter++) 
        {
            Event::McPartToTrajectoryRel* partToTraj = *trajIter;

            // Recover pointers to particle and trajectory
            const Event::McParticle*   part = partToTraj->getFirst();
            const Event::McTrajectory* traj = partToTraj->getSecond();

            // Get the particle properties...
            Event::McParticle::StdHepId hepid= part->particleProperty();
            ParticleProperty* ppty = m_ppsvc->findByStdHepID( hepid );
            if (ppty == 0) {
                log << MSG::DEBUG << "hepid = " << hepid << endreq;
                continue;
            }

            addTrack(traj->getPoints(), ppty->charge(), primary);
            primary=false;
        }
    } else {
        // use the endpoints to connect the dots
        SmartDataPtr<Event::McParticleCol>   mcPart(eventSvc(), "/Event/MC/McParticleCol");
        if(mcPart) for(Event::McParticleCol::const_iterator part=mcPart->begin()++; //note skip one
            part != mcPart->end(); part++) {

                Event::McParticle::StdHepId hepid= (*part)->particleProperty();
                ParticleProperty* ppty = m_ppsvc->findByStdHepID( hepid );
                if (ppty == 0) {
                  log << MSG::DEBUG << "hepid = " << hepid << endreq;
                  continue;
                }
                std::string name = ppty->particle(); 
                PointList line;
                Event::McTrajectoryPoint* initPoint = new Event::McTrajectoryPoint();
                initPoint->setPoint((*part)->initialPosition());
                line.push_back(initPoint);
                Event::McTrajectoryPoint* finalPoint = new Event::McTrajectoryPoint();
                finalPoint->setPoint((*part)->finalPosition());
                line.push_back(finalPoint);
                addTrack(line, ppty->charge(), primary);
                primary=false;

            }
    }
}




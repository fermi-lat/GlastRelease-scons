// $Header$
// Original Author: T. Burnett

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


// Detector stuff
#include "instrument/DetectorConverter.h"
#include "instrument/DetectorVolumes.h"
#include "instrument/Scintillator.h"
#include "instrument/CsIDetector.h"
#include "instrument/SiDetector.h"

// In order to set SiDetector parameters
#include "xml/IFile.h"

//gui
#include "GuiSvc/IGuiSvc.h"
#include "gui/DisplayControl.h"
#include "gui/GuiMgr.h"
#include "geomrep/BoxRep.h"


/*! \class IrfDisplay
\brief Display of data from the IRF

  */

class IrfDisplay : public Algorithm {

public:
  IrfDisplay(const std::string& name, ISvcLocator* pSvcLocator); 
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();

private:

    // the GlastDetSvc used for access to detector info
    IGlastDetSvc*    m_detSvc;
    
    // nested classes for display, that are also visitors
    class TileRep;
    TileRep* m_tileRep; // the ACD display rep
    class LogRep;
    LogRep*  m_logRep;
    class StripRep;
    StripRep* m_stripRep;
};

static const AlgFactory<IrfDisplay>  Factory;
const IAlgFactory& IrfDisplayFactory = Factory;
//#define HISTOGRAMS

#ifdef HISTOGRAMS
#include "GaudiKernel/IHistogram1D.h"
#include "Gaudikernel/IHistogramSvc.h"
static IHistogram1D  * h_strips, *h_tiles, *h_logs;
#endif

static bool m_shadow=false;  // use to turn detector outlines on/off

//------------------------------------------------------------------------------
IrfDisplay::IrfDisplay(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0) 
{
    declareProperty("detector_shadow",   m_shadow=false);
}

//------------------------------------------------------------------------------
//   Tile display Rep
//------------------------------------------------------------------------------
class   IrfDisplay::TileRep : public DetectorConverter , public gui::DisplayRep {
public:
    TileRep():m_count(0){}
    void   forward ( const Scintillator& d )
    {
        if ( d.empty()) return;
        if ( d.energy() < Scintillator::threshold()) return;

        const GlastDetector::_Volume * v = d.volume();
        const DetectorBox* db = dynamic_cast<const DetectorBox*> (d.volume());
        const Box& b = *db->box();
        setColor("red");
        append(BoxRep(b));
        ++m_count;
        setColor("black");
    }
    void clear(){
        m_count=0;
        DisplayRep::clear();
    }
    void update(){};
    int count()const { return m_count;}
private:
    int m_count;
};
//------------------------------------------------------------------------------
//   Log display rep
//------------------------------------------------------------------------------
class   IrfDisplay::LogRep : public DetectorConverter , public gui::DisplayRep {
public:
    LogRep():m_count(0){}
    double displayScale()const { return 100; }
    
    void   forward ( const CsIDetector& d )
    {
        if (d.empty()) return;
        const CoordTransform T(d.localToGlobal());  
        
        // First draw a shadow of the hit crystal.... 
        if(m_shadow) {
            const GlastDetector::_Volume * v = d.volume();
            const DetectorBox* db = dynamic_cast<const DetectorBox*> (d.volume());
            Box b(db->length(), db->width(), db->height()); // need to make a copy 
            b.transform(T);
            setColor("grey");  
            append(BoxRep(b));
        }

        float height = d.energy() * displayScale();
        if( height <= 0 ) return;
        
        // horizontal: make boxes with cross-section of the xtal,
        // lengths proportional to energy deposited. display at ends
        height = d.Lresp() * displayScale();
        if (height > 0.) {
            setColor("red");
            Box left_box(height,d.width(),d.depth());
            left_box.moveX(-(d.length()-height)/2);
            left_box.transform(T);
            append(BoxRep(left_box));
        }
        
        height = d.Rresp() * displayScale();
        if (height > 0.) {
            setColor("green");
            Box right_box(height,d.width(),d.depth());
            right_box.moveX((d.length()-height)/2);
            right_box.transform(T);
            append(BoxRep(right_box));
        }
        setColor("black");
        ++m_count;
    }
    void clear(){
        m_count=0;
        DisplayRep::clear();
    }
    void update(){ };
    int count()const { return m_count;}
private:
    int m_count;
};
//------------------------------------------------------------------------------
//   Strip display rep
//------------------------------------------------------------------------------
class   IrfDisplay::StripRep : public DetectorConverter , public gui::DisplayRep {
public:
    StripRep():m_count(0){}
    void   forward ( const SiDetector& d )
    {
        // static variable implmentation
        static float pScale = 0.175f;
        
        if ( d.empty()) return;
        CoordTransform T(d.localToGlobal());

        // First draw a shadow of the hit plane of silicon in grey
        if(m_shadow) {
            const GlastDetector::_Volume * v = d.volume();
            const DetectorBox* db = dynamic_cast<const DetectorBox*> (d.volume());

            Box b(db->length(), db->width(), db->height());// make a copy
            b.transform(T);
            setColor("grey");  
            append(BoxRep(b));
        }

        // Now draw the hits as red thermometers bars
        setColor("red");   
          
        for( SiDetector::const_iterator strip = d.begin(); strip != d.end(); ++strip) {
            float h = strip->energy() * pScale,
            x0 = d.localX(strip->index()),
            y0 = SiDetector::panel_width()/2+SiDetector::electronics_gap()/2;
            h *= 1.;
            Point from(x0,y0,0), to(x0,y0,h);
            from.transform(T); to.transform(T);
            moveTo(from);
            lineTo(to);
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
};

//------------------------------------------------------------------------------
/*! 
*/
StatusCode IrfDisplay::initialize() {
 using namespace gui;   
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();
    
    // now try to find the GlastDevSvc service
    StatusCode sc = service("GlastDetSvc", m_detSvc);

    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return sc;
    }

    // get the Gui service
    IGuiSvc* guiSvc=0;
    sc = service("GuiSvc", guiSvc);

    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GuiSvc!" << endreq;
        return sc;
    }

    // create new display and add the reps
    DisplayControl& display = guiSvc->guiMgr()->display();

    display.add(m_tileRep = new TileRep, "ACD hits");
    display.add(m_logRep = new LogRep,   "Log hits");
    display.add(m_stripRep = new StripRep, "Si Strip hits");
	display.menu().addSeparator();

    // SiDetector needs parameters from instrument/xml/instrument.xml to calculate strip positions
    SiDetector::loadParameters(* const_cast<xml::IFile*>(m_detSvc->iniFile()));
    // Scintillator needs parameters from instrument/xml/instrument.xml to determine veto threshold
    Scintillator::loadParams(* const_cast<xml::IFile*>(m_detSvc->iniFile()));
#ifdef HISTOGRAMS
        // book the histograms
    h_strips = histoSvc()->book("/stat/irf", 11, "StripCount", 50, 0, 50.);
    h_tiles  = histoSvc()->book("/stat/irf", 12, "TileCount", 50, 0, 50.);
    h_logs   = histoSvc()->book("/stat/irf", 13, "LogCount", 50, 0, 50.); 
#endif
    return sc;
}

//------------------------------------------------------------------------------
StatusCode IrfDisplay::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    m_detSvc->accept(*m_tileRep);
    m_detSvc->accept(*m_logRep);
    m_detSvc->accept(*m_stripRep);
    log << MSG::INFO << "Displaying " 
        << m_tileRep->count()  << " tiles, "
        << m_logRep->count()   << " logs, " 
        << m_stripRep->count() << " strips." << endreq;
#ifdef HISTOGRAMS
    h_strips->fill(m_stripRep->count(), 1.);
    h_logs->fill(m_logRep->count(), 1.);
    h_tiles->fill(m_tileRep->count(), 1.);
#endif
    return sc;
}


//------------------------------------------------------------------------------
StatusCode IrfDisplay::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}







/**
 * @file AcdDigiMergeAlg.cxx
 * @brief implementation  of the algorithm AcdDigiMergeAlg.
 *
 * @author Zach Fewtrell zachary.fewtrell@nrl.navy.mil
 * 
 *  $Header$
 */

// Gaudi specific include files
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

// Glast specific includes
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "AcdUtil/IAcdCalibSvc.h"

#include <map>

// Class definition
class AcdDigiMergeAlg : public Algorithm {

public:

    AcdDigiMergeAlg(const std::string& name, ISvcLocator* pSvcLocator); 

    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize() {return StatusCode::SUCCESS;}
 
private:

    /// used for constants & conversion routines.
    IGlastDetSvc*  m_detSvc;
  
    /// pointer to CalCalibSvc objects.
//    ICalCalibSvc*  m_overCalibSvc;  
    AcdUtil::IAcdCalibSvc*  m_digiCalibSvc;  
};


// Define the factory for this algorithm
static const AlgFactory<AcdDigiMergeAlg>  Factory;
const IAlgFactory& AcdDigiMergeAlgFactory = Factory;

/// construct object & declare jobOptions
AcdDigiMergeAlg::AcdDigiMergeAlg(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{

    // Declare the properties that may be set in the job options file
//    declareProperty("FirstRangeReadout",   m_firstRng= "autoRng"); 
}

/// initialize the algorithm. retrieve helper tools & services
StatusCode AcdDigiMergeAlg::initialize() 
{
    StatusCode sc = StatusCode::SUCCESS;
  
    MsgStream msglog(msgSvc(), name());
    msglog << MSG::INFO << "initialize" << endreq;

    sc = setProperties();
    if (sc.isFailure()) {
        msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
        return sc;
    }
    
    // try to find the GlastDetSvc service
    sc = service("GlastDetSvc", m_detSvc);
    if (sc.isFailure() ) {
        msglog << MSG::ERROR << "  can't get GlastDetSvc " << endreq;
        return sc;
    }
  
    // obtain the CalCalibSvc for the input overlay digis (data!)
//    sc = service("CalOverlayCalibSvc", m_overCalibSvc);
//    if (sc.isFailure()) 
//    {
//        msglog << MSG::ERROR << "can't get CalOverlayCalibSvc "  << endreq;
//        return sc;
//    }
  
    // now obtain the CalCalibSvc for the simulation digis
    sc = service("AcdSimCalibSvc", m_digiCalibSvc);
    if (sc.isFailure()) 
    {
        msglog << MSG::ERROR << "can't get CalCalibSvc "  << endreq;
        return sc;
    }

    return sc;
}

/// \brief take Hits from McIntegratingHits, create & register CalDigis
StatusCode AcdDigiMergeAlg::execute() 
{
    // Get a message object for output
    MsgStream msglog(msgSvc(), name());

    // First, recover any overlay digis, to see if we have anything to do
    SmartDataPtr<Event::AcdDigiCol> overlayDigiCol(eventSvc(), EventModel::Overlay::AcdDigiCol);
    if(!overlayDigiCol) return StatusCode::SUCCESS;

    // Now recover the digis for this event
    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(eventSvc(), EventModel::Digi::AcdDigiCol);

    // Create a map of the simulation digis, indexing by identifier
    std::map<idents::AcdId, Event::AcdDigi*> idToDigiMap;

    for(Event::AcdDigiCol::iterator acdIter = acdDigiCol->begin(); acdIter != acdDigiCol->end(); acdIter++)
    {
        Event::AcdDigi*     acdDigi = *acdIter;
        const idents::AcdId digiId  = acdDigi->getId();

        // Add this digi to our map
        idToDigiMap[digiId] = acdDigi;
    }

    // The calibration will need the Event Header information for both the sim and the overlay
    SmartDataPtr<Event::EventHeader> digiHeader(eventSvc(), "/Event");
  
    if (!digiHeader) {
      msglog << MSG::ERROR << "Unable to retrieve event timestamp for digis" << endreq;
      return StatusCode::FAILURE;
    }

    // The calibration will need the Event Header information for both the sim and the overlay
    SmartDataPtr<Event::EventHeader> overHeader(eventSvc(), EventModel::Overlay::EventHeader);
  
    if (!overHeader) {
      msglog << MSG::ERROR << "Unable to retrieve event timestamp for overlay" << endreq;
      return StatusCode::FAILURE;
    }

    // Loop through the input digis and using the above map merge with existing MC digis
    for(Event::AcdDigiCol::iterator overIter  = overlayDigiCol->begin(); overIter != overlayDigiCol->end(); overIter++)
    {
        Event::AcdDigi* overDigi = *overIter;

        const idents::AcdId overId = overDigi->getId();

        // We will need to correct the pedestals for the overlay digi, to match the simulation
        // The first step requires overwriting the time stamp information in the EventHeader, 
        // start by storing the current time locally
        TimeStamp currentTime(digiHeader->time());

        // Now overwrite with the time stamp from the overlay file
        digiHeader->setTime(overHeader->time());

        // Return the time stamp information to normal
        digiHeader->setTime(currentTime);

        // Does this overlay digi match one in the sim map?
        std::map<idents::AcdId,Event::AcdDigi*>::iterator acdIter = idToDigiMap.find(overId);

        // If this digi is not already in the CalDigiCol then no merging needed, just add it
        if (acdIter == idToDigiMap.end())
        {
            Event::AcdDigi* digi = new Event::AcdDigi();
            *digi = *overDigi;
            digi->setParent(0);
//            digi->addToStatus(Event::CalDigi::DIGI_OVERLAY);

            // Ok do something, anything, here

            acdDigiCol->push_back(digi);
        }
        // Otherwise, we need to merge the digi
        else
        {
            // This would be where one did something to combine two digis... ugh
            int j = 0;
        }
    }

    return StatusCode::SUCCESS;
}

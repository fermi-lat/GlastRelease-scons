#include "CalDigiAlg.h"
#include "CalUtil/CalFailureModeSvc.h"
/// Gaudi specific include files
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IToolSvc.h"

/// Glast specific includes
#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "CalibData/CalibTime.h"

// Relational Table
#include "Event/RelTable/Relation.h"
#include "Event/RelTable/RelTable.h"
#include "Event/Digi/CalDigi.h"
#include "CLHEP/Random/RandGauss.h"
/// for min and floor functions
#include <algorithm>
#include <math.h>
#include <iostream.h>

// to access an XML containing Digi parameters file
#include "xml/IFile.h"

// std stuff
#include <utility>
#include <string>

// Define the factory for this algorithm
static const AlgFactory<CalDigiAlg>  Factory;
const IAlgFactory& CalDigiAlgFactory = Factory;



// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.


CalDigiAlg::CalDigiAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    // Declare the properties that may be set in the job options file
    declareProperty ("xmlFile", m_xmlFile="$(CALDIGIROOT)/xml/CalDigi.xml");
    declareProperty ("taperToolName", m_taperToolName="LinearTaper");
    declareProperty ("convertAdcToolName", 
                                m_convertAdcToolName="LinearConvertAdc");
    declareProperty ("doFluctuations", m_doFluctuations="yes");
    declareProperty( "startTime", m_startTimeAsc = "2003-1-10_00:20");
    declareProperty( "calibFlavor", m_calibFlavor = "none");
}

StatusCode CalDigiAlg::initialize() {
    // Purpose and Method: initialize the algorithm. Set up parameters from detModel
    // Inputs: detModel parameters
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
  
    // extracting int constants from detModel. Store into local cache - member variables.
   
    double value;
    typedef std::map<int*,std::string> PARAMAP;

    PARAMAP param;
    param[&m_xNum]=        std::string("xNum");
    param[&m_yNum]=        std::string("yNum");
    param[&m_eTowerCal]=   std::string("eTowerCAL");
    param[&m_eLatTowers]=  std::string("eLATTowers");
    param[&m_CalNLayer]=   std::string("CALnLayer");
    param[&m_nCsIPerLayer]=std::string("nCsIPerLayer");
    param[&m_nCsISeg]= std::string("nCsISeg");
    param[&m_eXtal]=       std::string("eXtal");
    param[&m_eDiodeMSmall]=std::string("eDiodeMSmall");
    param[&m_eDiodePSmall]=std::string("eDiodePSmall");
    param[&m_eDiodeMLarge]=std::string("eDiodeMLarge");
    param[&m_eDiodePLarge]=std::string("eDiodePLarge");
    param[&m_eMeasureX]=std::string("eMeasureX");
    param[&m_eMeasureY]=std::string("eMeasureY");
    param[m_noise]=std::string("cal.noiseLarge");
    param[m_noise+1]=std::string("cal.noiseSmall");
    param[m_ePerMeV+1]=std::string("cal.ePerMeVSmall");
    param[m_ePerMeV]=std::string("cal.ePerMevLarge");
    param[&m_pedestal]=std::string("cal.pedestal");
    param[&m_maxAdc]=std::string("cal.maxAdcValue");
   

    // now try to find the GlastDevSvc service
   
    IGlastDetSvc* detSvc;
    StatusCode sc = service("GlastDetSvc", detSvc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to get GlastDetSvc " << endreq;
        return sc;
    }
    
    for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
        if(!detSvc->getNumericConstByName((*it).second, &value)) {
            log << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
            return StatusCode::FAILURE;
        } else *((*it).first)=(int)value;
    }
    
   
    // extracting double constants from detModel. Store into local cache - member variables.
    
    typedef std::map<double*,std::string> DPARAMAP;
    DPARAMAP dparam;
    dparam[m_maxEnergy]=std::string("cal.maxResponse0");
    dparam[m_maxEnergy+1]=std::string("cal.maxResponse1");
    dparam[m_maxEnergy+2]=std::string("cal.maxResponse2");
    dparam[m_maxEnergy+3]=std::string("cal.maxResponse3");
    dparam[&m_lightAtt]=std::string("cal.lightAtt");
    dparam[&m_CsILength]=std::string("CsILength");
    dparam[&m_thresh]=std::string("cal.zeroSuppressEnergy");
    
    for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
        if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
            log << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
            return StatusCode::FAILURE;
        } 
    }
    
    
    // Read in the parameters from the XML file
    xml::IFile m_ifile(m_xmlFile.c_str());
    if (m_ifile.contains("cal","ePerMeVinDiode")) 
        m_ePerMeVinDiode = m_ifile.getDouble("cal", "ePerMeVinDiode");
    else return StatusCode::FAILURE;
  
  
    sc = toolSvc()->retrieveTool(m_taperToolName,m_taper);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_taperToolName << endreq;
        return sc;
    }
   
    sc = toolSvc()->retrieveTool(m_convertAdcToolName,m_convertAdc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to create " << m_convertAdcToolName << endreq;
        return sc;
    }

    m_doFluctuationsBool = (m_doFluctuations == "yes") ? 1 : 0;
   
    sc = service("CalFailureModeSvc", m_FailSvc);
    if (sc.isFailure() ) {
        log << MSG::INFO << "  Did not find CalFailureMode service" << endreq;
  m_FailSvc = 0;
    }
   
    // get pointer to CalibDataSvc
    sc = service("CalibDataSvc", m_pCalibDataSvc, true);

    if ( !sc.isSuccess() ) {
      log << MSG::ERROR 
      << "Could not get IDataProviderSvc interface of CalibDataSvc" 
      << endreq;
      return sc;
    }
 
    // Query the IDetDataSvc interface of the calib data service
    sc = m_pCalibDataSvc->queryInterface(IID_IDetDataSvc, 
                                (void**) &m_detDataSvc);
    if ( !sc.isSuccess() ) {
      log << MSG::ERROR 
      << "Could not query IDetDataSvc interface of CalibDataSvc" 
      << endreq;
      return sc;
    } else {
      log << MSG::DEBUG 
      << "Retrieved IDetDataSvc interface of CalibDataSvc" 
      << endreq;
    }
   
    // Get properties from the JobOptionsSvc
    sc = setProperties();
    if ( !sc.isSuccess() ) {
      log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
      return sc;
     }

    unsigned int underpos = m_startTimeAsc.find("_");
    if (underpos < m_startTimeAsc.size()) {
         m_startTimeAsc.replace(underpos, 1, " ");
    }
    m_startTime = facilities::Timestamp(m_startTimeAsc).getClibTime();

    log << MSG::DEBUG << "Properties were read from jobOptions" << endreq;
    log << MSG::INFO << "Time of first event: (ascii) "
        << m_startTimeAsc       << endreq; 
    log << MSG::INFO << "Time of first event: (seconds since 1970) "
        << m_startTime       << endreq; 
   
    return StatusCode::SUCCESS;
}


StatusCode CalDigiAlg::execute() {

    // Purpose and Method: take Hits from McIntegratingHit and perform the following steps:
    //  for deposit in a crystal segment, take into account light propagation to the two ends and 
    // apply light taper based on position along the length.
    //  keep track of direct deposit in the diode.
    // Then combine diode (with appropriate scale factor) and crystal deposits and add noise.
    // Then, add noise to 'unhit' crystals.
    //  Finally, convert to ADC units and pick the range for hits above threshold.
    // Inputs: McIntegratingHit
    
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    
    //Take care of insuring that data area has been created
    DataObject* pNode = 0;
    sc = eventSvc()->retrieveObject( EventModel::Digi::Event /*"/Event/Digi"*/, pNode);
    
    if (sc.isFailure()) {
        sc = eventSvc()->registerObject(EventModel::Digi::Event /*"/Event/Digi"*/,new Event::DigiEvent);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register " << EventModel::Digi::Event /*<< /Event/Digi "*/ << endreq;
            return sc;
        }
    }
    
    // Set the event time
    facilities::Timestamp time = facilities::Timestamp(m_startTime);
    //log << MSG::DEBUG << "Event time: "
        //<< time.getString()
        //<< endreq; 
    CalibData::CalibTime ctime(time);
    //log << MSG::DEBUG << "Event time (hours) " << ctime.hours() << endreq;
    m_detDataSvc->setEventTime(ctime);

    DataObject *pObject;

    // setup calibration information
    pPeds = 0;
    pGains = 0;

    if(m_calibFlavor != "none") {
      std::string fullPedPath = "/Calib/CAL_Ped/"+m_calibFlavor;
      std::string fullGainPath = "/Calib/CAL_ElecGain/"+m_calibFlavor;

      //getting pointers to the calibration data of each type
      if(m_pCalibDataSvc->retrieveObject(fullPedPath, pObject) 
          == StatusCode::SUCCESS) {
        pPeds = dynamic_cast<CalibData::CalCalibPed *> (pObject);
        if (!pPeds) {
          log << MSG::ERROR << "Dynamic cast to CalCalibPed failed" << endreq;
          return StatusCode::FAILURE;
        }
      } else
        log << MSG::INFO 
            << "Unable to retrieve pedestals from calib database" << endreq;

      if(m_pCalibDataSvc->retrieveObject(fullGainPath, pObject) 
          == StatusCode::SUCCESS) {
        pGains = dynamic_cast<CalibData::CalCalibGain *> (pObject);
        if (!pGains) {
          log << MSG::ERROR << "Dynamic cast to CalCalibGain failed" << endreq;
          return StatusCode::FAILURE;
        }
      } else            
        log << MSG::INFO 
            << "Unable to retrieve gains from calib database" << endreq;
    } 

    //  clear signal array: map relating xtal signal to id. Map holds diode and crystal responses
    //  separately during accumulation.
    
    m_signalMap.clear();
    m_idMcInt.clear();
    
    sc = fillSignalEnergies();
    if (sc != StatusCode::SUCCESS) return sc;
    
    //sc = addNoiseToSignals();
    //if (sc != StatusCode::SUCCESS) return sc;
    
    //sc = addNewNoiseHits();
    //if (sc != StatusCode::SUCCESS) return sc;
    
    
    
    log << MSG::DEBUG; if (log.isActive()){ log.stream() << m_signalMap.size() << "calorimeter hits in m_signalMap";}  log << endreq;
    for( SignalMap::iterator jit=m_signalMap.begin(); jit!=m_signalMap.end();jit++){
        log << MSG::DEBUG; if (log.isActive()){ log.stream() << " id " << (*jit).first
            << " s0=" << (*jit).second.getSignal(idents::CalXtalId::POS)
            << " s1=" << (*jit).second.getSignal(idents::CalXtalId::NEG)
            << " d0=" << (*jit).second.getDiodeEnergy(idents::CalXtalId::POS)
            << " d1=" << (*jit).second.getDiodeEnergy(idents::CalXtalId::NEG)
            << " d2=" << (*jit).second.getDiodeEnergy(idents::CalXtalId::POS+2)
            << " d3=" << (*jit).second.getDiodeEnergy(idents::CalXtalId::NEG+2)
           ;} log << endreq;
    }
    
    sc = createDigis();
    if (sc != StatusCode::SUCCESS) return sc;    
    
    return sc;
}

StatusCode CalDigiAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}

StatusCode CalDigiAlg::createDigis() {
    // Purpose and Method: 
    // create digis on the TDS from the deposited energies

    MsgStream log(msgSvc(), name());

    Event::CalDigiCol* digiCol = new Event::CalDigiCol;
    
    Event::RelTable<Event::CalDigi, Event::McIntegratingHit> digiHit;
    digiHit.init();
    
    // iterate through the SignalMap to look at ALL Xtals and have CalUtil add 
    // electronic noise
    // if either side is above threshold, then select the appropriate ADC range and
    // create a readout
    for (int tower = 0; tower < m_xNum*m_yNum; tower++)
      for (int layer = 0; layer < m_CalNLayer; layer++)
        for (int col = 0; col < m_nCsIPerLayer; col++){
          char rangeP,rangeM;
          unsigned short adcP, adcM;
          unsigned short status = 0;
          idents::CalXtalId xtalId(tower,layer,col);
           
          // loop over the plus and minus faces of the crystal 
          double resp[8];
          for (int face=0; face<2; face++){
            unsigned short adc;
            int range;
    
            // resp is the signal in the 2 diodes on that face.
            SignalMap::iterator nit=  m_signalMap.find(xtalId);
            if( nit != m_signalMap.end() ){
              XtalSignal* signal = &(*nit).second;
              resp[0+4*face] = signal->getDiodeEnergy(face
                  +2*idents::CalXtalId::LARGE);
              resp[2+4*face] = signal->getDiodeEnergy(face
                  +2*idents::CalXtalId::SMALL);
              resp[1+4*face] = resp[0+4*face];
              resp[3+4*face] = resp[2+4*face];
            } else 
              resp[0+4*face]=resp[2+4*face]=resp[1+4*face]=resp[3+4*face]=0;

            // check for failure mode. If killed, set to zero and set DEAD bit
            if ((m_FailSvc != 0) &&  
              (m_FailSvc->matchChannel(xtalId,
              (idents::CalXtalId::XtalFace)face))) {

              range = 0;
              status = (face == idents::CalXtalId::POS) ?  
                (status | Event::CalDigi::CalXtalReadout::DEAD_P) : 
                (status | Event::CalDigi::CalXtalReadout::DEAD_N);

            } else {
              // calculate the ADC values, adding the noise, from the deposited
              // energy in the two diodes on this face and find the best range
              range= m_convertAdc->calculateAdcAndNoise(xtalId,
                  (idents::CalXtalId::XtalFace)face, 
                  resp+4*face, pPeds, pGains );

              // resp[range]==0 if and only if LEX8<zero_suppress
              // in which case there is no hit on this end of the log
              if( resp[range+4*face]==0 ) continue;

              adc= (short unsigned int) resp[range+4*face];
            }

            // assign the plus/minus readouts
            if(face == idents::CalXtalId::POS){ adcP=adc; rangeP=range;}
            else { adcM=adc; rangeM=range;}
          }

          log << MSG::DEBUG; 
          if (log.isActive()){ 
            log.stream() <<" id=" << xtalId 
            << " rangeP=" << int(rangeP) << " adcP=" << adcP
            << " rangeM=" << int(rangeM) << " adcM=" << adcM;
          } 
          log << endreq;
           

          // set status to ok for POS and NEG if no other bits set.
      
          if ((status & 0x00FF) == 0) status = 
              (status | Event::CalDigi::CalXtalReadout::OK_P);
          if ((status & 0xFF00) == 0) status = 
              (status | Event::CalDigi::CalXtalReadout::OK_N);
        
          // set up the digi
          Event::CalDigi::CalXtalReadout read = 
            Event::CalDigi::CalXtalReadout(rangeP, adcP, rangeM, adcM, status);
          Event::CalDigi* curDigi = 
            new Event::CalDigi(idents::CalXtalId::BESTRANGE, xtalId);
          curDigi->addReadout(read);
          //Event::CalDigi* curDigi= new Event::CalDigi(
                                    //idents::CalXtalId::ALLRANGE, xtalId);
          //for( int range =0; range<4; range ++ ){
            //adcP= (short int) resp[range];
            //adcM= (short int) resp[range+4];
            //Event::CalDigi::CalXtalReadout read = 
              //Event::CalDigi::CalXtalReadout( range, adcP, range, adcM, status);
            //curDigi->addReadout(read);
          //}
          
          // set up the relational table between McIntegratingHit and digis
          typedef std::multimap< idents::CalXtalId, Event::McIntegratingHit* >::const_iterator ItHit;
          std::pair<ItHit,ItHit> itpair = m_idMcInt.equal_range(xtalId);
            
          for (ItHit mcit = itpair.first; mcit!=itpair.second; mcit++)
          {
            Event::Relation<Event::CalDigi,Event::McIntegratingHit> *rel =
            new Event::Relation<Event::CalDigi,Event::McIntegratingHit>(curDigi,mcit->second);
            digiHit.addRelation(rel);
          }
            
          // add the digi to the digi collection
          digiCol->push_back(curDigi);
        }
    
    StatusCode sc = eventSvc()->registerObject(EventModel::Digi::CalDigiCol, digiCol);
    
    if (!(sc == StatusCode::FAILURE))
        sc = eventSvc()->registerObject(EventModel::Digi::CalDigiHitTab,digiHit.getAllRelations());
    return sc;
}

StatusCode CalDigiAlg::fillSignalEnergies() {
    // Purpose and Method: collect deposited energies from McIntegratingHits and store
    // in map sorted by XtalID. 
    // multimap used to associate mcIntegratingHit to id. There can be multiple
    // hits for the same id.  
    
    // get McIntegratingHit collection. Abort if empty.
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    SmartDataPtr<Event::McIntegratingHitVector> McCalHits(eventSvc(),EventModel::MC::McIntegratingHitCol ); //"/Event/MC/IntegratingHitsCol");
    
    if (McCalHits == 0) {
        log << MSG::DEBUG; if (log.isActive()){ log.stream() << "no calorimeter hits found" ;} log << endreq;
        return sc;
    }
    
    // loop over hits - pick out CAL hits
    
    for (Event::McIntegratingHitVector::const_iterator it = McCalHits->begin(); it != McCalHits->end(); it++) {
        
        
        //   extracting hit parameters - get energy and first moment
        
        idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*it)->volumeID());
        
        
        //   extracting parameters from volume Id identifying as in CAL
        
        if ((int)volId[fLATObjects] == m_eLatTowers &&
            (int)volId[fTowerObjects] == m_eTowerCal){ 
            
            double ene = (*it)->totalEnergy();
            HepPoint3D mom1 = (*it)->moment1();
            
            //  log << MSG::DEBUG <<  "McIntegratingHits info \n"  
            //  << " ID " << volId.name()
            //  <<  " energy " << ene
            //  <<  " moments " << mom1.x()
            //  << endreq;
            
            int col = volId[fCALXtal];
            int layer = volId[fLayer];
            int towy = volId[fTowerY];
            int towx = volId[fTowerX];
            int tower = m_xNum*towy+towx; 
            int segm = volId[fSegment];
            
            idents::CalXtalId mapId(tower,layer,col);
            
            log << MSG::DEBUG; if (log.isActive()){ log.stream() <<  "Identifier decomposition \n"  
                << " col " << col
                << " layer " << layer
                << " towy " << towy
                << " towx " << towx
                << " segm " << segm
                ;} log << endreq;
            
            XtalSignal& xtalSignalRef = m_signalMap[mapId];
            
            // Insertion of the id - McIntegratingHit pair
            
            m_idMcInt.insert(std::make_pair(mapId,*it));
            
            if((int)volId[fCellCmp] == m_eDiodePLarge)
                xtalSignalRef.addDiodeEnergy(ene,idents::CalXtalId::POS);
            
            else if((int)volId[fCellCmp] == m_eDiodeMLarge)
                xtalSignalRef.addDiodeEnergy(ene,idents::CalXtalId::NEG);
            
            else if((int)volId[fCellCmp] == m_eDiodePSmall )
                xtalSignalRef.addDiodeEnergy(ene,idents::CalXtalId::POS+2);
            
            else if((int)volId[fCellCmp] ==  m_eDiodeMSmall) 
                xtalSignalRef.addDiodeEnergy(ene,idents::CalXtalId::NEG+2);
            
            else if((int)volId[fCellCmp] ==  m_eXtal ){
                
                
                // let's define the position of the segment along the crystal
                
                double relpos = (segm+0.5)/m_nCsISeg;
                
                
                // in local reference system x is always oriented along the crystal
                double dpos =  mom1.x(); 
                
                relpos += dpos/m_CsILength;
                
                
                
                // take into account light tapering
                
                std::pair<double,double> signals;
                signals = m_taper->calculateSignals( mapId, relpos, ene );
                
                // set up a XtalMap to act as key for the map - if there is no entry, add
                // add one, otherwise add signal to existing map element.
                
                
                xtalSignalRef.addSignal(signals.first,signals.second);
                
                
                
            }
        }
    }   
    return sc;
}

StatusCode CalDigiAlg::addNoiseToSignals() {
    // Purpose and Method: 
    // add electronic noise to the diode response and add to the crystal 
    // response. The diodeEnergy becomes the readout source.
    
    for(SignalMap::iterator mit=m_signalMap.begin(); mit!=m_signalMap.end();mit++){
        XtalSignal* xtal_signal = &(*mit).second;
        
        for ( int idiode =0; idiode < 4; idiode++){
            int diode_type = idiode/2;
            int face  = idiode%2;
            double signal = xtal_signal->getSignal(face);
            
            double diode_ene = xtal_signal->getDiodeEnergy(idiode);
            
            
            // convert energy deposition in a diode to
            // the equivalent energy in a crystal 
            diode_ene *= m_ePerMeVinDiode/m_ePerMeV[diode_type];
            
            // add crystal signal - now diode energy contains
            // the signal at given diode in energy units
            // (equivalent energy deposition at the crystal center)
            diode_ene += signal;
            
            // add poissonic fluctuations in the number of electrons in a diode
            if (m_doFluctuationsBool) {
                float numberElectrons = diode_ene * m_ePerMeV[diode_type];
                
                // approximate Poisson distribution by gaussian 
                // for numberElectrons >> 1
                float electronFluctuation = sqrt(numberElectrons) 
                  * RandGauss::shoot();
                
                diode_ene += electronFluctuation /m_ePerMeV[diode_type];
            }
            
             //add electronic noise: now done in CalUtil
             //diode_ene += RandGauss::shoot()
               //*m_noise[diode_type]/m_ePerMeV[diode_type];
            
            //store modified diode signal in the signal map
            xtal_signal->setDiodeEnergy(diode_ene,idiode);
            
        }
    }
    return StatusCode::SUCCESS;
}

/*
StatusCode CalDigiAlg::addNewNoiseHits() {
    // Purpose and Method: 
    //  adding electronic noise to the channels without signal,
    //   storing it in the signal map if  one of crystal faces
    //  has the noise above the threshold   

    // noise in MeV for the large diode
    double noise_MeV = double(m_noise[idents::CalXtalId::LARGE])
                      /double(m_ePerMeV[idents::CalXtalId::LARGE]); 
    for (int tower = 0; tower < m_xNum*m_yNum; tower++){
        for (int layer = 0; layer < m_CalNLayer; layer++){
            for (int col = 0; col < m_nCsIPerLayer; col++){

                double eneM = noise_MeV*RandGauss::shoot();
                double eneP = noise_MeV*RandGauss::shoot();
                idents::CalXtalId mapId(tower,layer,col);

                if((eneM > m_thresh || eneP > m_thresh) &&
                    m_signalMap.find(mapId) == m_signalMap.end()){
                    XtalSignal& xtalSignalRef = m_signalMap[mapId];
                    xtalSignalRef.addDiodeEnergy(eneM,idents::CalXtalId::POS);
                    xtalSignalRef.addDiodeEnergy(eneP,idents::CalXtalId::NEG);

                }
            }
        }
    }
    return StatusCode::SUCCESS;;
}*/

CalDigiAlg::XtalSignal::XtalSignal() {
    // Purpose and Method: default constructor setting signals to zero
    m_signal[idents::CalXtalId::POS] = 0;
    m_signal[idents::CalXtalId::NEG] = 0;
    
    for(int i=0; i<4; i++)m_Diodes_Energy.push_back(0.);
    
}

CalDigiAlg::XtalSignal::XtalSignal(double s1, double s2) {
    // Purpose and Method: constructor setting signals to those input
    m_signal[idents::CalXtalId::POS] = s1;
    m_signal[idents::CalXtalId::NEG] = s2;
    
    for(int i=0; i<4; i++)m_Diodes_Energy.push_back(0.);
    
}
void CalDigiAlg::XtalSignal::addSignal(double s1, double s2) {
    // Purpose and Method: add signals s1, s2 to already existing signals.
    m_signal[idents::CalXtalId::POS] += s1;
    m_signal[idents::CalXtalId::NEG] += s2;
    return;
} 

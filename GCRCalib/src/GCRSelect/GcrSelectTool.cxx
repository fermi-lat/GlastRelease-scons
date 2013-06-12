#include "IGcrSelectTool.h"
#include "DataStructures.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartRefVector.h"

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


#include "TkrUtil/ITkrGeometrySvc.h"

#include "enums/TriggerBits.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/ThreeVector.h" 
#include "CLHEP/Vector/LorentzVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/GcrReconClasses.h"
#include "Event/Recon/CalRecon/GcrSelectClasses.h"
#include "Event/TopLevel/Event.h"

#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McParticle.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"


#include "geometry/Vector.h"
#include "geometry/Ray.h"

#include "TMath.h"
#include "TObjArray.h"
#include "TObject.h"

#include "OnboardFilterTds/ObfFilterStatus.h"

/**   
 * @class GcrSelectTool
 *
 */

//-----------------------------------------------------------------------------------------------------------------
class GcrSelectTool : public IGcrSelectTool,  public AlgTool 
{
    public:    
        GcrSelectTool(const std::string & type, const std::string & name, const IInterface * parent );
        virtual ~GcrSelectTool() {};

        /// @brief Intialization of the tool
        virtual StatusCode initialize();
        /// @brief Default cluster finding framework
        virtual StatusCode GcrSelectTool::selectGcrXtals();


    private:
        /// PRIVATE METHODS

        DataSvc* getDataSvc(){return m_dataSvc;} 
        void setDataSvc(DataSvc* dataSvc){m_dataSvc = dataSvc;}

        ///loads Cal dimensions parameters: 
        StatusCode readGlastDet();


        Event::GcrXtalVec getGcrXtalVec(){return m_gcrXtalVec;} 
        void setGcrXtalVec(Event::GcrXtalVec gcrXtalVec){m_gcrXtalVec = gcrXtalVec;}
        Event::GcrReconVals*   m_gcrReconVals; 



        float GcrSelectTool::chooseEth(); // decides energy threshold to select useful hits for calibration, depending of which OBF filter (Gamma, HFC, Mip, DFC) the event passes


        StatusCode retrieveGcrXtalsMap();

        void attributeGcrGrade ();

        void findClusters();

        void fillExpectedEdep();
        void verifyExpectedEdep();

        int inferZ(double correctedEnergy);

        int buildHitsMap();
        void verifyHitsMap();


        int buildSelectedXtalsVec();

        double GcrSelectTool::correctEnergy(idents::CalXtalId xtalId, Event::GcrXtal gcrXtal);

        StatusCode storeGcrSelectedXtals ();
        StatusCode storeGcrSelectVals () ;

        //StatusCode buildOutputs();
        StatusCode store();


        //The ideal thing would be to write a method findClusters() that would be a mix of buildCluArr and fillDataVectors
        void buildEnergyMap(); // this method should not really be necessary, as information is already contained in m_hitsMap

        void buildLayMultArr();

        void buildCluArr();
        void displayCluArr();

        void fillDataVectors();

        void verifyDataVectors();

        void buildLayMatchTrackCrit();
        void verifyLayMatchTrackCrit();
        void buildLayMatchMultCrit();
        void verifyLayMatchMultCrit();


        /// PRIVATE DATA MEMBERS

        static const int NTOW = 16;
        static const int NLAY = 8;
        static const int NCOL = 12;

        /// MsgStream member variable to speed up execution
        MsgStream          m_log;

        /// Pointer to the Gaudi data provider service
        DataSvc*           m_dataSvc;

        /// Pointer to the Geant4 propagator
        IPropagator * m_G4PropTool; 

        /// the GlastDetSvc used for access to detector info
        IGlastDetSvc*      m_detSvc;

        /// TkrGeometrySvc used for access to tracker geometry info
        ITkrGeometrySvc*   m_geoSvc;

        //Xtal dimensions:
        double             m_xtalHeight,m_xtalWidth,m_xtalLength;

        //data structures:


        /// collection of GcrXtals. Each GcrXtal contains a pointer
        /// to corresponding CalXtalRecData, if any, and the value
        /// of the pathlength defined by first MCParticle direction in
        /// the Xtal. Pointer points
        /// to null if no CalXtalRecData is found.
        Event::GcrXtalVec m_gcrXtalVec;     


        /// collection of GcrSelectedXtals. Each GcrXtal contains a pointer
        /// to corresponding CalXtalRecData, raw-energy and corrected energy values, 
        /// selection grade, and the value
        /// of the pathlength defined by first MCParticle direction in
        /// the Xtal.
        Event::GcrSelectedXtalsVec m_gcrSelectedXtalsVec;     

        /// collection of GcrXtals, used to store info into TDS
        Event::GcrSelectedXtalsCol*   m_gcrSelectedXtalsCol; 
        Event::GcrSelectVals*   m_gcrSelectVals; 


        unsigned int m_gcrOBFStatusWord;  // contains information about Gamma, HFC, Mip, DFC filters vetoes
        //{'Gam':0,'Hfc':1,'Mip':2,'Dfc':3}


        ZVect m_zVect; 

        //data structures necessary for E.N. codes:
        /// Hits (Xtals with energy deposit) Map, pointer points on null if no XtalRecData for Xtal  
        Event::CalXtalRecData*  m_hitsMap[NTOW][NLAY][NCOL]; 

        /// Xtals Map, contains -1 if no energy deposit is expected from MCparticle trajectory extrapolation
        int  m_gcrXtalsMap[NTOW][NLAY][NCOL]; 

        GcrTowersVec m_gcrTowersVec;

        int m_layerMultiplicity[NTOW][NLAY];  //multiplicity per layer
        float m_cluArr[NTOW][NLAY][NCOL][4]; //third dimension is NCOL because there is a max of NCOL clusters per layer
        //[O]:total energy in cluster
        //[1]:first Xtal index in cluster
        //[2]:last Xtal index in cluster
        //[3]:number of hits in cluster

        double m_energyMap[NTOW][NLAY][NCOL];

        float m_Eth; //Energy threshold for considering an energy deposit as representatif

        int m_inferedZ; /// inferedZ for this event, from cluster from first lay.

        // global Criteria match variables: intertowers criteria
        int m_layMatchTrackCrit[8];
        int m_layMatchMultCrit[8];
        int m_layMatchMultCrit2[8];
        //allows display if debugging:
        bool m_debugging;



} ;


//-----------------------------------------------------------------------------------------------------------------
//static ToolFactory<GcrSelectTool> s_factory;
//const IToolFactory& GcrSelectToolFactory = s_factory;
DECLARE_TOOL_FACTORY(GcrSelectTool);

//-----------------------------------------------------------------------------------------------------------------
GcrSelectTool::GcrSelectTool(const std::string & type, 
        const std::string & nameIn,
        const IInterface * parent ) : AlgTool( type, nameIn, parent ),
    m_log(msgSvc(), name())
{ 
    declareInterface<IGcrSelectTool>(this) ; 


    return;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode GcrSelectTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    m_log.setLevel(outputLevel());
    m_debugging=false;

    m_log << MSG::INFO << "GcrSelectTool BEGIN initialize()" << endreq ;


    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    // find TkrGeometrySvc service
    if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
        m_log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // find G4 propagation tool
    if(!toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
        m_log << MSG::ERROR << "Couldn't find the G4PropationTool!" << endreq;
        return StatusCode::FAILURE;
    }

    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
        m_log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }


    m_log << MSG::INFO << "GcrSelectTool END initialize()" << endreq ;  

    ///loads Cal dimensions parameters: 
    readGlastDet();


    return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode GcrSelectTool::readGlastDet()
{
    StatusCode sc = StatusCode::SUCCESS;
    //m_log << MSG::INFO << "GcrSelectTool BEGIN readGlastDet()" << endreq ;  


    //Xtal dimensions:
    m_detSvc->getNumericConstByName("CsIHeight", &m_xtalHeight);  
    m_detSvc->getNumericConstByName("CsIWidth", &m_xtalWidth);  
    m_detSvc->getNumericConstByName("CsILength", &m_xtalLength);  
    //m_log << MSG::INFO << "CsIHeight,CsIWidth,CsILength=" << m_xtalHeight << "," << m_xtalWidth<< "," << m_xtalLength<< endreq ;  


    //m_log << MSG::INFO << "GcrSelectTool END readGlastDet()" << endreq ;  
    return sc;
}



//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------

StatusCode GcrSelectTool::selectGcrXtals(){
    m_log << MSG::INFO << "GcrSelectTool BEGIN selectGCRXtals()" << endreq ; 

    StatusCode sc = StatusCode::SUCCESS; 

    if (retrieveGcrXtalsMap() == StatusCode::SUCCESS){



        buildEnergyMap();

        findClusters();

        fillExpectedEdep();


        fillDataVectors(); 

        buildLayMatchTrackCrit();

        if(m_debugging)
            verifyLayMatchTrackCrit();

        buildLayMatchMultCrit();

        if(m_debugging)
            verifyLayMatchMultCrit();


        // for DEBUG: verifyExpectedEdep();


        if(m_debugging)
            verifyDataVectors();


        // buildOutputs();

        if(buildSelectedXtalsVec() >0)
            store();

    }

    m_log << MSG::INFO << "GcrSelectTool END selectGCRXtals()" << endreq ;  

    return sc; 


}

//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::buildLayMatchTrackCrit(){

    /// projection of m_gcrXtalsMap on LAYERS
    bool debugging = false;

    if(debugging)
        m_log << MSG::INFO << "BEGIN buildLayMatchTrackCrit in GcrSelectTool" << endreq;

    for(int itow=0; itow<NTOW;itow++)  
        for (int ilay=0;ilay<NLAY;ilay++)
            for (int icol=0;icol<NCOL;icol++)
                m_layMatchTrackCrit[ilay] = -1;


    for(int ilay=0; ilay<NLAY;ilay++) { 

        if(debugging)
            m_log << MSG::INFO << "ilay= " << ilay << endreq;

        for (int itow=0;itow<NTOW;itow++)
        {
            for (int icol=0;icol<NCOL;icol++)
            {
                //m_log << MSG::INFO << "m_gcrXtalsMap[" << itow << "][" << ilay << "][" << icol << "]= " << m_gcrXtalsMap[itow][ilay][icol]<< endreq;
                //m_log << MSG::INFO << "m_hitsMap[" << itow << "][" << ilay << "][" << icol << "]= NULL? " << (m_hitsMap[itow][ilay][icol] == NULL) << endreq;
                if(m_gcrXtalsMap[itow][ilay][icol] != -1){  // if the track crosses this Xtal
                    m_layMatchTrackCrit[ilay] = 1;
                    if(debugging)
                        m_log << MSG::INFO << "m_gcrXtalsMap[" << itow << "][" << ilay << "][" << icol << "]= " << m_gcrXtalsMap[itow][ilay][icol]<< endreq;
                }	       
            }// end of loop on columns
        }// end of loop on towers
    }//end of loop on layers


    // tag as not good Layers, all layers beyond the first not good layer 
    bool fstGoodLayFound = false;
    for(int ilay=0; ilay<NLAY;ilay++) { 

        if(m_layMatchTrackCrit[ilay]==1){
            if(!fstGoodLayFound) fstGoodLayFound = true;
        }
        else{
            if(fstGoodLayFound){
                for(int i=ilay; i<NLAY; i++)
                    m_layMatchTrackCrit[i]=-2;

                return;
            }

        }
    }   


    if(debugging)
        m_log << MSG::INFO << "END buildLayMatchTrackCrit in GcrSelectTool" << endreq;
}

void GcrSelectTool::verifyLayMatchTrackCrit(){

    m_log << MSG::INFO << "BEGIN verifyLayMatchTrackCrit in GcrSelectTool" << endreq;
    //display:

    for (int ilay=0;ilay<NLAY;ilay++)
        m_log << MSG::INFO << "m_layMatchTrackCrit[" << ilay << "] = " << m_layMatchTrackCrit[ilay] << endreq;

    m_log << MSG::INFO << "END verifyLayMatchTrackCrit in GcrSelectTool" << endreq;


}


//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::buildLayMatchMultCrit(){


    m_log << MSG::VERBOSE << "BEGIN buildLayMatchMultCrit in GcrSelectTool" << endreq;

    for(int itow=0; itow<NTOW;itow++)  
        for (int ilay=0;ilay<NLAY;ilay++)
            for (int icol=0;icol<NCOL;icol++){
                m_layMatchMultCrit[ilay] = -1;
                m_layMatchMultCrit2[ilay] = -1;
            }

    for(int ilay=0; ilay<NLAY;ilay++) { 
        int nbClu=0;
        int maxNbHits=0;
        m_log << MSG::VERBOSE << "ilay= " << ilay << endreq;

        for (int itow=0;itow<NTOW;itow++)
        {
            bool lastArrayFound =false;

            for(int iClu=0; iClu<NCOL; iClu++){

                int nbHitsInCluster=0; 
                if(m_cluArr[itow][ilay][iClu][0]>0)
                {
                    m_log << MSG::VERBOSE << "m_cluArr[" << itow<<"]["<<ilay<<"]["<<iClu<<"][0]=" << m_cluArr[itow][ilay][iClu][0] << endreq;
                    nbClu++;

                    m_log << MSG::VERBOSE << "m_cluArr[" << itow<<"]["<<ilay<<"]["<<iClu<<"][1]," << "m_cluArr[" << itow<<"]["<<ilay<<"]["<< iClu <<"][2] = " << m_cluArr[itow][ilay][iClu][1]<<","<< m_cluArr[itow][ilay][iClu][2]<< endreq;	     

                    nbHitsInCluster =  (int)(m_cluArr[itow][ilay][iClu][2] -  m_cluArr[itow][ilay][iClu][1] + 1);

                    m_log << MSG::VERBOSE << "nbHitsInCluster= " << nbHitsInCluster << endreq;



                    if(nbHitsInCluster > maxNbHits)
                        maxNbHits = nbHitsInCluster;

                }
                else
                    lastArrayFound=true;

            }


        }// end of loop on towers

        m_log << MSG::VERBOSE << "nbClu= " << nbClu << endreq;
        m_log << MSG::VERBOSE << "maxNbHits= " << maxNbHits << endreq;

        if(nbClu>0)
        {
            if(nbClu<2 && maxNbHits<3)
                m_layMatchMultCrit[ilay] = 1;
            else if(nbClu>=2 && maxNbHits>=3)
                m_layMatchMultCrit[ilay] = 4;
            else if(nbClu>=2)
                m_layMatchMultCrit[ilay] = 2;
            else if(maxNbHits>=3)
                m_layMatchMultCrit[ilay] = 3;
        }


    }//end of loop on layers

    // tag as not good Layers, all layers beyond the first not good layer 
    bool fstLayWithClusters = false;
    for(int ilay=0; ilay<NLAY;ilay++) { 

        if(m_layMatchMultCrit[ilay]!=-1)// if clusters have been found on this layer
            if(!fstLayWithClusters) fstLayWithClusters = true;

        if(m_layMatchMultCrit[ilay]==1){ // if layer is good
            m_layMatchMultCrit2[ilay]=m_layMatchMultCrit[ilay];
            m_log << MSG::VERBOSE << "m_layMatchMultCrit[ " << ilay << "]==1"<< endreq;

        }
        else{   // if layer is not good
            m_log << MSG::VERBOSE << "m_layMatchMultCrit[ " << ilay << "]<>1"<< endreq;

            int layMatchMultCritCourant = m_layMatchMultCrit[ilay];
            if(m_layMatchMultCrit[ilay] == -1 && !fstLayWithClusters)  // layer is not good but no cluster has been found before 
                continue;

            else{// layer is not good, clusters have been found, every layer beyond this one has to be tagged as bad one
                for(int i=ilay; i<NLAY; i++)
                    m_layMatchMultCrit2[i]=layMatchMultCritCourant;

                return;
            }    

        }

    }   

    m_log << MSG::VERBOSE << "END buildLayMatchMultCrit in GcrSelectTool" << endreq;
}

void GcrSelectTool::verifyLayMatchMultCrit(){

    m_log << MSG::VERBOSE << "BEGIN verifyLayMatchMultCrit in GcrSelectTool" << endreq;
    //display:

    for (int ilay=0;ilay<NLAY;ilay++)
        m_log << MSG::DEBUG << "m_layMatchMultCrit[" << ilay << "] = " << m_layMatchMultCrit[ilay] << endreq;

    for (int ilay=0;ilay<NLAY;ilay++)
        m_log << MSG::DEBUG << "m_layMatchMultCrit2[" << ilay << "] = " << m_layMatchMultCrit2[ilay] << endreq;

    m_log << MSG::DEBUG << "END verifyLayMatchMultCrit in GcrSelectTool" << endreq;


}



//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::buildEnergyMap()
{ 

    bool debugging = false;
    if(debugging)
        m_log << MSG::INFO << "BEGIN buildEnergyMap in GcrSelectTool" << endreq;

    int numHits=0;
    int itow, ilay, icol;

    //initialization of Xtal(Hits) Ids 
    for (itow=0; itow<16; itow++)
        for (ilay=0; ilay<8; ilay++)
            for (icol=0; icol<12; icol++){
                m_energyMap[itow][ilay][icol]=-1000.0;		
            }



    Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 


    if (calXtalRecCol){
        // loop over CalXtalRecdata
        //m_log << MSG::DEBUG << "juste before loop, calXtalRecCol->numberOfObjects()= " << calXtalRecCol->numberOfObjects() << endreq;
        for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
        {

            Event::CalXtalRecData* xTalData = *xTalIter;

            itow=xTalData->getPackedId().getTower();
            ilay=xTalData->getPackedId().getLayer();
            icol=xTalData->getPackedId().getColumn();

            /**Event::CalMipTrack& calMipTrack = *trackIter;
            // Need to create a new CalMipTrack which will be "owned" by the TDS
            Event::CalMipTrack* newMipTrack = new Event::CalMipTrack();

            // Now copy to this new track. This should properly copy all elements (including inherited vector)
             *newMipTrack = calMipTrack;*/


            m_energyMap[itow][ilay][icol] = xTalData->getEnergy();

            if(debugging){
                m_log << MSG::INFO << "CalXtalRecData found at:"<< itow<<"/"<<ilay<<"/"<<icol<<endreq; 
                //m_log << MSG::INFO << "Energy="<< m_energyMap[itow][ilay][icol] <<endreq;
            } 


            //if(debugging)
            //m_log << MSG::INFO << "m_energyMap["<< itow<<"]["<<ilay<<"]["<<icol<<"]= "<< 
            //m_log << MSG::INFO << "m_hitsMap[itow][ilay][icol]->getPosition()= " << m_hitsMap[itow][ilay][icol]->getPosition() << endreq;
            //(m_hitsMap[itow][ilay][icol])->getPosition().printOn(std::cout);
            // m_log << MSG::INFO << "\n" << endreq;

            numHits++;

        }
    }
    /**else 
      m_log << MSG::INFO << "no calXtalRecCol " << endreq;
     */	
    //DEBUG: display
    if(debugging){
        for (itow=0; itow<16; itow++)
            for (ilay=0; ilay<8; ilay++)
                for (icol=0; icol<12; icol++){
                    if(m_energyMap[itow][ilay][icol]>0)
                        m_log << MSG::INFO << "m_energyMap[" << itow << "][" << ilay << "][" << icol <<"]= " << m_energyMap[itow][ilay][icol] << endreq;
                }
    }

    if(debugging)
        m_log << MSG::INFO << "END buildEnergyMap in GcrSelectTool, numHits= " << numHits << endreq;
    //return numHits;
}

//------------------------------------------------------------------------------------------------------------


StatusCode GcrSelectTool::retrieveGcrXtalsMap(){
    //m_log << MSG::INFO << "GcrSelectTool BEGIN retrieveMcGcrXtalMap()" << endreq ; 

    StatusCode sc = StatusCode::SUCCESS;

    /// collection of GcrXtals, used to store info into TDS
    Event::GcrXtalCol*   gcrXtalCol; 



    int numGcrXtals =0; 

    int itow = -1;
    int ilay = -1;
    int icol = -1;

    //initialization of gcrXtalsMap
    for (int iTow=0; iTow<NTOW; iTow++)
        for (int iLay=0; iLay<NLAY; iLay++)
            for (int iCol=0; iCol<NCOL; iCol++)
                m_gcrXtalsMap[iTow][iLay][iCol]=-1;

    m_gcrXtalVec.clear();

    //retrieve gcrXtalCol
    gcrXtalCol = SmartDataPtr<Event::GcrXtalCol>(m_dataSvc, EventModel::CalRecon::GcrXtalCol);

    //update of gcrXtalsMap & gcrXtalVec
    if(gcrXtalCol){ 
        //m_log << MSG::INFO << "McGcrXtalCol found" << endreq ; 
        for(Event::GcrXtalCol::const_iterator gcrXtalColIter=gcrXtalCol->begin(); gcrXtalColIter != gcrXtalCol->end(); gcrXtalColIter++)
        {
            Event::GcrXtal* currentGcrXtal= *gcrXtalColIter;
            //Event::CalXtalRecData* xTalData = currentGcrXtal->getXtal(); // this is not NULL for sure, as no gcrXtals are stored by GcrRecon
            // when no corresponding XtalRecData is found
            idents::CalXtalId xtalId = currentGcrXtal->getXtalId();

            itow=xtalId.getTower();
            ilay=xtalId.getLayer();
            icol=xtalId.getColumn();

            //m_log << MSG::INFO << "xTalData found: (" << itow <<","<< ilay <<","<< icol << ")= (" << itow <<"," << ilay << "," << icol <<")" << endreq;
            m_gcrXtalsMap[itow][ilay][icol] = numGcrXtals++;

            m_gcrXtalVec.push_back(Event::GcrXtal(currentGcrXtal->getXtalId(), currentGcrXtal->getPathLength(), currentGcrXtal->getClosestFaceDist(), currentGcrXtal->getCrossedFaces(), currentGcrXtal->getEntryPoint(), currentGcrXtal->getExitPoint()));
            //m_log << MSG::INFO << "m_gcrXtalVec.size()= " << m_gcrXtalVec.size() << endreq ;

        }

        //DEBUG: display McGcrXtalsMap
        /**for (int itow=0; itow<NTOW; itow++)
          for (int ilay=0; ilay<NLAY; ilay++)
          for (int icol=0; icol<NCOL; icol++)
          if(m_gcrXtalsMap[itow][ilay][icol]>=0)
          m_log << MSG::INFO << "m_gcrXtalsMap[" << itow <<"]["<< ilay <<"]["<< icol << "]= " << m_gcrXtalsMap[itow][ilay][icol]<< endreq;
         */
    }
    else{
        sc = StatusCode::FAILURE;
        m_log << MSG::INFO << "no GcrXtalCol found" << endreq ; 
    }


    //m_log << MSG::INFO << "GcrSelectTool END retrieveMcGcrXtalMap()" << endreq ;  
    return sc;

}


//-----------------------------------------------------------------------------------------------------------------
StatusCode GcrSelectTool::store(){

    StatusCode sc = StatusCode::SUCCESS;


    storeGcrSelectedXtals();
    storeGcrSelectVals();

    return sc;
}


//-----------------------------------------------------------------------------------------------------------------


void GcrSelectTool::attributeGcrGrade (){
    // m_log << MSG::INFO << "GcrSelectTool BEGIN attributeGCRGrade()" << endreq ;  
    // m_log << MSG::INFO << "GcrSelectTool END attributeGCRGrade()" << endreq ;  


}

//-----------------------------------------------------------------------------------------------------------------


void GcrSelectTool::findClusters(){

    if(m_debugging)
        m_log << MSG::INFO << "GcrSelectTool BEGIN findClusters()" << endreq ;


    buildLayMultArr();

    buildCluArr();  

    //DEBUG: 
    /**if(m_debugging)
      displayCluArr();*/

    if(m_debugging)
        m_log << MSG::INFO << "GcrSelectTool END findClusters()" << endreq ;  
}


//-----------------------------------------------------------------------------------------------------------------


float GcrSelectTool::chooseEth(){

    float Eth=100.0;
    bool passGamma=0, passHIP=0, passMIP=0, passDGN=0;

    //gcrOBFStatusWord comes from gcrReconVals - calculated in GcrReconTool

    Event::GcrReconVals* gcrReconVals;
    gcrReconVals = SmartDataPtr<Event::GcrReconVals>(m_dataSvc,EventModel::CalRecon::GcrReconVals);
    if(!(gcrReconVals == 0)){
        m_gcrOBFStatusWord = gcrReconVals->getGcrOBFStatusWord();
    }else{
        m_log << MSG::INFO << "No gcrReconVals found" << endreq;
    }

    //   passGamma = (m_gcrOBFStatusWord & (1 << OnboardFilterTds::ObfFilterStatus::GammaFilter))!=0;
    //   passHIP = (m_gcrOBFStatusWord   & (1 << OnboardFilterTds::ObfFilterStatus::HIPFilter))!=0;
    //   passMIP = (m_gcrOBFStatusWord   & (1 << OnboardFilterTds::ObfFilterStatus::MIPFilter))!=0;
    //   passDGN = (m_gcrOBFStatusWord   & (1 << OnboardFilterTds::ObfFilterStatus::DGNFilter))!=0;

    passGamma = m_gcrOBFStatusWord & 1;
    passHIP = m_gcrOBFStatusWord   & 2;
    passMIP = m_gcrOBFStatusWord   & 4;
    passDGN = m_gcrOBFStatusWord   & 8;

    m_log << MSG::DEBUG << "passGamma, passHIP, passMIP, passDGN= " << passGamma <<"," << passHIP << "," <<passMIP << "," << passDGN << endreq;

    SmartDataPtr<Event::EventHeader> pEvent(m_dataSvc, EventModel::EventHeader);
    unsigned int word2 =  ( pEvent==0? 0 : pEvent->triggerWordTwo());
    unsigned int Trig_gemengine = ((word2 >> enums::ENGINE_offset) & enums::ENGINE_mask);
    bool engine4ON = (Trig_gemengine==4);

    if (passHIP)
        Eth = 100.0;
    else if (passGamma) {
        if(engine4ON)
        { 
            Eth = 100.0;
        } 
        else 
        {
            Eth = 5.0;
        }
    }
    else if (passMIP || passDGN)
    {
        Eth = 5.0;
    }

    m_log << MSG::DEBUG << "Eth =" << Eth << endreq; 

    return Eth;

}


//-----------------------------------------------------------------------------------------------------------------


void GcrSelectTool::buildLayMultArr(){

    m_log << MSG::VERBOSE << "GcrSelectTool BEGIN buildLayMultArr()" << endreq ; 

    m_Eth = chooseEth();  //Energy threshold for considering an energy deposit as representatif



    for (int itow=0; itow<NTOW; itow++){
        //m_log << MSG::INFO << "itow=" << itow << endreq;

        for (int ilay=0; ilay<NLAY; ilay++){

            //m_log << MSG::INFO << "ilay=" << ilay << endreq;

            int multCurrentLay=0;
            for (int icol=0; icol<NCOL; icol++){
                // m_log << MSG::INFO << "m_energyMap[itow][ilay][icol]" << m_energyMap[itow][ilay][icol] << endreq;
                if(m_energyMap[itow][ilay][icol]>=m_Eth) {
                    //m_log << MSG::INFO << "multiplicity increment" << endreq;
                    multCurrentLay++;

                }
            }//end of for icol=0...
            m_layerMultiplicity[itow][ilay] = multCurrentLay;
        }//end of for int ilay=0...
    }//end of for int itow=0...

    //DISPLAY:
    /**if(m_debugging){
      for (int itow=0; itow<NTOW; itow++)
      for (int ilay=0; ilay<NLAY; ilay++)
      if(m_layerMultiplicity[itow][ilay]>0)
      m_log << MSG::INFO << "m_layerMultiplicity[" << itow<< "][" << ilay<< "]=" << m_layerMultiplicity[itow][ilay] << endreq;
      }*/

    m_log << MSG::VERBOSE << "GcrSelectTool END buildLayMultArr()" << endreq ; 

}


//-----------------------------------------------------------------------------------------------------------------


void GcrSelectTool::buildCluArr(){
    //m_log << MSG::INFO << "GcrSelectTool BEGIN buildCluArr()" << endreq ; 


    for(int itow=0; itow<NTOW;itow++) { //Car ici on a toutes les tours.

        for (int ilay=0;ilay<NLAY;ilay++)
        {

            for (int nclu=0;nclu<NCOL;nclu++)
            {
                m_cluArr[itow][ilay][nclu][0]= 0.; //total energy in cluster
                m_cluArr[itow][ilay][nclu][1]=-1;  //first xtal index in cluster
                m_cluArr[itow][ilay][nclu][2]=-1;  //last  xtal index in cluster
                m_cluArr[itow][ilay][nclu][3]=-1;  //number of hits in cluster

            }



            int nhl=m_layerMultiplicity[itow][ilay];//nb of hits in this layer

            if(nhl<1) //no hit in this layer
                continue;
            else
            {
                int nha=0;
                int nhb=0;

                int nclu=0;
                while((nha<NCOL)&&(nhb<NCOL))
                {


                    while( (m_energyMap[itow][ilay][nha]<m_Eth) && (nha<NCOL))
                    {
                        //m_log << MSG::INFO << "Dans while energy<m_Eth, nha= " << nha  << ", m_energyMap[itow][ilay][nha]= " << m_energyMap[itow][ilay][nha] << ", seuilE=" << seuilE << endreq;
                        nha++;

                    }

                    //m_log << MSG::INFO << "à la sortie de while En<m_Eth nha= " << nha<< endreq;
                    if(nha>=NCOL)
                        continue;  

                    //m_log << MSG::INFO << "itow,ilay,nclu=" << itow << "," << ilay << "," << nclu << endreq;

                    nclu++; //at least one cluster found

                    if(m_cluArr[itow][ilay][nclu][1]<0)
                        m_cluArr[itow][ilay][nclu][1]=nha; //first xtal index in this cluster

                    //look for the last hit of the cluster:
                    nhb=nha;

                    while( (m_energyMap[itow][ilay][nhb]>=m_Eth) && (nhb<NCOL) )  
                    {
                        m_cluArr[itow][ilay][nclu][0]+=m_energyMap[itow][ilay][nhb];	
                        nhb++;

                    }

                    //m_log << MSG::INFO << "à la sortie de while En>=m_Eth nhb= " << nhb<< endreq;
                    if(nhb > NCOL)
                        continue;

                    //m_log << MSG::INFO << "m_cluArr[itow][ilay][nclu][2]= " << m_cluArr[itow][ilay][nclu][2]<< endreq;
                    if(m_cluArr[itow][ilay][nclu][2]<0)
                    {
                        //m_log << MSG::INFO << "dans m_cluArr[itow][ilay][nclu][2]<0"<< endreq;
                        m_cluArr[itow][ilay][nclu][2]=nhb-1; //last xtal index in this cluster
                        m_cluArr[itow][ilay][nclu][3]=nhb-nha; //number of hits in this cluster
                        //(attention en sortie de boucle on a nhb++)
                    }                      


                    //nclu++;
                    nha=nhb;
                }//end of while(nha<NCOL && nhb<NCOL)
            }//end of else (hit found for this layer)
        }//end of ilay=0...

    }//end of itow=0....

    //m_log << MSG::INFO << "GcrSelectTool END buildCluArr()" << endreq ; 

}


//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::displayCluArr()
{ 
    //m_log << MSG::INFO << "GcrSelectTool BEGIN displayCluArr()" << endreq ; 
    //DISPLAY RESULTS
    for (int itow=0; itow<NTOW; itow++)
        for (int ilay=0; ilay<NLAY; ilay++)
            for (int iclu=0; iclu<NCOL; iclu++){
                if(m_cluArr[itow][ilay][iclu][0] > 0.){
                    m_log << MSG::INFO << "m_cluArr[" << itow<< "][" << ilay<< "][" << iclu <<"][0]=" << m_cluArr[itow][ilay][iclu][0] << endreq;
                    m_log << MSG::INFO << "m_cluArr[" << itow<< "][" << ilay<< "][" << iclu <<"][1]=" << m_cluArr[itow][ilay][iclu][1] << endreq;
                    m_log << MSG::INFO << "m_cluArr[" << itow<< "][" << ilay<< "][" << iclu <<"][2]=" << m_cluArr[itow][ilay][iclu][2] << endreq;
                    m_log << MSG::INFO << "m_cluArr[" << itow<< "][" << ilay<< "][" << iclu <<"][3]=" << m_cluArr[itow][ilay][iclu][3] << endreq;
                }
            }

    //m_log << MSG::INFO << "GcrSelectTool END displayCluArr()" << endreq ; 
}

//-----------------------------------------------------------------------------------------------------------------
/**
 *  This method builts a 3D map[NTOW][NLAY][NCOL], whose entries are pointers to XTalRecData, if any
 */
int GcrSelectTool::buildHitsMap()
{ 

    //m_log << MSG::INFO << "BEGIN buildHitsMap in GcrSelectTool" << endreq;

    int numHits=0;
    int itow, ilay, icol;

    //initialization of Xtal(Hits) Ids 
    for (itow=0; itow<NTOW; itow++)
        for (ilay=0; ilay<NLAY; ilay++)
            for (icol=0; icol<NCOL; icol++){
                m_hitsMap[itow][ilay][icol]=0;		
            }



    Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 


    if (calXtalRecCol){
        // loop over CalXtalRecdata
        //m_log << MSG::DEBUG << "juste before loop, calXtalRecCol->numberOfObjects()= " << calXtalRecCol->numberOfObjects() << endreq;
        for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
        {

            Event::CalXtalRecData* xTalData = *xTalIter;

            itow=xTalData->getPackedId().getTower();
            ilay=xTalData->getPackedId().getLayer();
            icol=xTalData->getPackedId().getColumn();

            m_hitsMap[itow][ilay][icol] = xTalData;

            // m_log << MSG::INFO << "m_hitsMap[itow][ilay][icol]->getPosition()= " << endreq;
            //(m_hitsMap[itow][ilay][icol])->getPosition().printOn(std::cout);
            // m_log << MSG::INFO << "\n" << endreq;

            numHits++;

        }
    }
    else 
        m_log << MSG::INFO << "no calXtalRecCol " << endreq;



    //m_log << MSG::INFO << "END buildHitsMap in GcrSelectTool, numHits= " << numHits << endreq;
    return numHits;
}

// ----------------------------------------------------------------------------
void GcrSelectTool::verifyHitsMap(){

    //m_log << MSG::INFO << "BEGIN verifyHitsMap in GcrSelectTool"<< endreq;

    for (int itow=0; itow<NTOW; itow++)
        for (int ilay=0; ilay<NLAY; ilay++)
            for (int icol=0; icol<NCOL; icol++){
                Event::CalXtalRecData* xTalData = m_hitsMap[itow][ilay][icol];
                if(xTalData)
                    m_log << MSG::INFO << "m_hitsMap[" << itow << "][" << ilay << "][" << icol << "]->getPackedId()= " << xTalData->getPackedId()<< endreq;

            }

    //m_log << MSG::INFO << "END verifyHitsMap in GcrSelectTool"<< endreq;
}




//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::fillDataVectors(){

    if(m_debugging)
        m_log << MSG::INFO << "GcrSelectTool BEGIN fillDataVectors()" << endreq ;  


    buildHitsMap();
    //verifyHitsMap();

    m_gcrTowersVec.clear();

    //########creation of GcrTowers and GcrLayers vector entries:
    for (int itow=0; itow<NTOW; itow++){//loop on towers

        //create towers vector entry:
        m_gcrTowersVec.push_back(GcrTower(itow));

        // create layers vector:
        GcrLayersVec gcrLayersVec;
        for (int ilay=0; ilay<NLAY; ilay++){//loop on layers
            //create layers vector entry:
            gcrLayersVec.push_back(GcrLayer(ilay,0,false,false));
        }

        //, and put it into gcrTowersVec
        GcrTower& currentGcrTower = m_gcrTowersVec.back();
        currentGcrTower.setLayersVec(gcrLayersVec);

    }// end of loop on towers


    //########creation of GcrClustersVectors and association to layers, only if clusters found.
    for (int itow=0; itow<NTOW; itow++){//loop on towers
        for (int ilay=0; ilay<NLAY; ilay++){//loop on layers

            //Variables linked to filtering criteria:
            // a cluster satisfies the clusterPredTrack criterium if it contains at most 2 hits and there is AT LEAST ONE of them on the predicted track
            // a layer satisfies the layerPredTrack criterium if EVERY cluster associated to it satisfies predTrack

            bool isGoodCluster; // does the cluster satisfies predictedTrackCriterium?
            bool isGoodLayer = true;// does the layer satisfies predictedTrackCriterium?

            GcrLayersVec& p_gcrLayerVec = m_gcrTowersVec.at(itow).getLayersVec();

            GcrClustersVec* gcrClustersVec=NULL;

            for (int iclu=0; iclu<NCOL; iclu++){//loop on clusters

                isGoodCluster = false; 
                //m_log << MSG::INFO << "iclu= " << iclu << endreq;

                if(m_cluArr[itow][ilay][iclu][0]>0.){  // if total energy in cluster >0...

                    // we create or pickup a gcrClustersVec
                    if ( (p_gcrLayerVec.at(ilay).getClustersVec()) &&
!(p_gcrLayerVec.at(ilay).getClustersVec()->empty()) ) // if a ClustersVec has already been associated to this layer..
                        gcrClustersVec = p_gcrLayerVec.at(ilay).getClustersVec(); // we take this clustersVec already present
                    else
                        gcrClustersVec = new GcrClustersVec();// else, we create a new one.



                    double totalCluCorrEnergy = 0.0;

                    int firstCluHitNb = (int)m_cluArr[itow][ilay][iclu][1];
                    int lastCluHitNb = (int)m_cluArr[itow][ilay][iclu][2];

                    GcrHitsVec gcrHitsVec;

                    //m_log << MSG::INFO << "itow,ilay,iclu= " << itow << "," << ilay << "," << iclu << "| firstCluHitNb,lastCluHitNb = " << firstCluHitNb << "," << lastCluHitNb << endreq ;

                    //build gcrHitsVec
                    double xtalCorrEnergy;


                    //int inferedZ=0;
                    for(int i=firstCluHitNb; i<=lastCluHitNb; i++){
                        xtalCorrEnergy= -999.0;
                        int mapXtal = m_gcrXtalsMap[itow][ilay][i];// if!=-1, there is a corresponding GcrXtal [itow][ilay][i]	    

                        bool isGoodHit = (mapXtal != -1); // true if hit is on the predicted track

                        if(m_debugging)
                            m_log << MSG::INFO << "itow,ilay,i: correspondingGcrXtalExist? " << itow << "," << ilay << ","<< i << ":" << (mapXtal!=-1) << endreq;


                        if(mapXtal != -1){
                            Event::GcrXtal currentGcrXtal= m_gcrXtalVec.at(mapXtal);
                            xtalCorrEnergy = correctEnergy(currentGcrXtal.getXtalId(), currentGcrXtal);

                            //INFER Z
                            /**if((xtalCorrEnergy>0) && (ilay==0))
                              inferedZ = inferZ(xtalCorrEnergy); /// ATTENTION: Il ne faut pas calculer le Z par Xtal, mais par cluster.
                            /// Voir diagrammes des traces sur le papier 

                             */
                            totalCluCorrEnergy+=xtalCorrEnergy;
                        }

                        isGoodCluster = isGoodCluster? isGoodCluster:isGoodHit;


                        gcrHitsVec.push_back(GcrHit(m_hitsMap[itow][ilay][i], isGoodHit, xtalCorrEnergy));


                    }// end of loop on currentClusters hits


                    isGoodCluster = (gcrHitsVec.size() <=2)? isGoodCluster:false;// if the cluster contains more than two hits, this criterium is not satisfied

                    gcrClustersVec->push_back(GcrCluster(gcrHitsVec,totalCluCorrEnergy, isGoodCluster));//we calculate total cluster corrected energy from individual hits corrected energy

                    p_gcrLayerVec.at(ilay).setClustersVec(gcrClustersVec);


                    isGoodLayer = isGoodLayer? isGoodCluster:isGoodLayer;// if false, it must remain false, if true, it depends on if currentCluster matches criterium or not
                }//end of if(m_cluArr[][][])

            }//end of loop on clusters

            if(!gcrClustersVec)
                isGoodLayer=false;
            else if(gcrClustersVec->size()!=1)
                isGoodLayer=false;

            //m_log << MSG::INFO << "isGoodLayer?" << isGoodLayer << ", gcrClustersVec exists?" << (gcrClustersVec!=NULL) << endreq ; 
            //if(gcrClustersVec)  m_log << MSG::INFO << "gcrClustersVec->size()=" << gcrClustersVec->size() << endreq ; 

            p_gcrLayerVec.at(ilay).setIsGoodLayer(isGoodLayer);


        }//end of loop on layers

    }//end of loop on towers

    if(m_debugging) 
        m_log << MSG::INFO << "GcrSelectTool END fillDataVectors()" << endreq ;  
}


//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::verifyDataVectors (){
    m_log << MSG::INFO << "GcrSelectTool BEGIN verifyDataVectors()" << endreq ;  

    //m_log << MSG::INFO << "m_gcrTowersVec.size()= " << m_gcrTowersVec.size() << endreq;
    int itow=0;
    for(GcrTowersVec::iterator gcrTowersVecIter = m_gcrTowersVec.begin(); gcrTowersVecIter != m_gcrTowersVec.end(); gcrTowersVecIter++){
        GcrTower currentTow = *gcrTowersVecIter;
        GcrLayersVec& currentLayersVec = currentTow.getLayersVec();
        //m_log << MSG::INFO << "itow:"<< itow << " gcrLayersVec.size()= " << currentLayersVec->size() << endreq;

        int ilay=0;
        for(GcrLayersVec::iterator gcrLayersVecIter= currentLayersVec.begin();
gcrLayersVecIter!= currentLayersVec.end(); gcrLayersVecIter++){

            GcrLayer currentLay = *gcrLayersVecIter;
            m_log << MSG::INFO << "itow,ilay:"<< itow << "," << ilay << " currentLay.getIsGoodLayer() " << currentLay.getIsGoodLayer() << " " << endreq;

            GcrClustersVec* currentClustersVec = currentLay.getClustersVec();
            if(currentClustersVec){
                m_log << MSG::INFO << "itow,ilay=" << itow << "," << ilay << " currentClustersVec.size()= " << currentClustersVec->size() << endreq;

                int iclu=0;

                for(GcrClustersVec::iterator currentClustersVecIter = currentClustersVec->begin(); currentClustersVecIter != currentClustersVec->end(); currentClustersVecIter++){
                    GcrCluster currentCluster = *currentClustersVecIter;

                    GcrHitsVec& currentHitsVec = currentCluster.getHitsVec();

                    if(currentHitsVec.size()>0){
                        m_log << MSG::INFO << "itow,ilay,iclu=" << itow << "," << ilay << "," << iclu << " " << endreq;
                        m_log << MSG::INFO << "itow,ilay,iclu=" << itow << ","
<< ilay << "," << iclu << " gcrHitsVec.size()= " << currentHitsVec.size() << endreq;
                        m_log << MSG::INFO << "cluster Total Corr Energy= " << currentCluster.getTotalCorrectedEnergy() << endreq;
                        //m_log << MSG::INFO << "ilay:"<< ilay << " currentLay.getIsGoodLayer() " << currentLay.getIsGoodLayer() << " "<< endreq;
                        m_log << MSG::INFO << "isGoodCluster= " << currentCluster.getIsGoodCluster() << " " << endreq;

                        int currentHitNb = 0;
                        for(GcrHitsVec::iterator currentHitsVecIter =
currentHitsVec.begin(); currentHitsVecIter != currentHitsVec.end(); currentHitsVecIter++){
                            GcrHit currentHit = *currentHitsVecIter;
                            //m_log << MSG::INFO << "XtalpackedId, iHit=" << currentHit.getXtalData()->getPackedId()<< ","<< currentHitNb << " currentHit corr energy= " << currentHit.getCorrectedEnergy() << endreq;
                            m_log << MSG::INFO << "XtalpackedId, iHit=" << currentHit.getXtalData()->getPackedId()<< ","<< currentHitNb << " currentHit isGoodHit= " << currentHit.getIsGoodHit() << " " << endreq;
                            //<< "inferedZ="<< currentHit.getInferedZ()<<endreq;
                            currentHitNb++;
                        }
                    }else //impossible case normally
                        m_log << MSG::INFO << "no hitsVec associated to itow,ilay,iclu=" << itow << "," << ilay << "," << iclu  << endreq;


                    iclu++;

                }// end ofloop on ClustersVec

            }// end of if(currentClustersVec)
            /**else{
              m_log << MSG::INFO << "no cluster vec associated to itow,ilay=" << itow << "," << ilay << endreq;

              }*/
            ilay++;

        }// end of loop on LayersVec
        itow++;

    }//end of loop on TowersVec


    m_log << MSG::INFO << "GcrSelectTool END verifyDataVectors()" << endreq ;  
}


//-----------------------------------------------------------------------------------------------------------------

double GcrSelectTool::correctEnergy(idents::CalXtalId xtalId, Event::GcrXtal gcrXtal){

    // m_log << MSG::INFO << "GcrSelectTool BEGIN correctEnergy()" << endreq ; 

    double correctedEnergy = 999.999;

    // find CalXtalRecData Energy
    //double depositedEnergy = xTalData -> getEnergy();
    double depositedEnergy = m_energyMap[xtalId.getTower()][xtalId.getLayer()][xtalId.getColumn()];
    // find gcrXtal pathLength
    double pathLength = gcrXtal.getPathLength();

    //m_log << MSG::INFO << "pathLength= " << pathLength << " ,xtalHeight= " << xtalHeight << endreq ; 

    //Calculate correctedEnergy;
    correctedEnergy = (depositedEnergy / pathLength) * m_xtalHeight;

    //m_log << MSG::INFO << "correctedEnergy= " << correctedEnergy << endreq ; 

    // m_log << MSG::INFO << "GcrSelectTool END correctEnergy()" << endreq ; 

    return correctedEnergy;




}

//-----------------------------------------------------------------------------------------------------------------


void GcrSelectTool::fillExpectedEdep(){
    //m_log << MSG::INFO << "GcrSelectTool BEGIN fillExpectedEdep()" << endreq ;  

    m_zVect.clear();

    //double minE, maxE;
    for(int z = 3; z<=26; z++){// loop on atomic numbers
        SpEnergiesInterval* spEnergiesInterval = new SpEnergiesInterval();
        if(z==3){

            spEnergiesInterval->setMinE(100.0);
            spEnergiesInterval->setMaxE(150.0);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }else     if(z==4){

            spEnergiesInterval->setMinE(150.0);
            spEnergiesInterval->setMaxE(270.0);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }else     if(z==5){

            spEnergiesInterval->setMinE(270.0);
            spEnergiesInterval->setMaxE(350.0);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }else     if(z==6){

            spEnergiesInterval->setMinE(350.0);
            spEnergiesInterval->setMaxE(527.333);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }else if(z==7){

            spEnergiesInterval->setMinE(527.333);
            spEnergiesInterval->setMaxE(703.198);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }

        else if(z==8){

            spEnergiesInterval->setMinE(703.198);
            spEnergiesInterval->setMaxE(898.573);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }

        else if(z==9){

            spEnergiesInterval->setMinE(898.573);
            spEnergiesInterval->setMaxE(1112.31);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }

        else if(z==10){

            spEnergiesInterval->setMinE(1112.31);
            spEnergiesInterval->setMaxE(1368.84);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==11){

            spEnergiesInterval->setMinE(1368.84);
            spEnergiesInterval->setMaxE(1627.78);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==12){

            spEnergiesInterval->setMinE(1627.78);
            spEnergiesInterval->setMaxE(1934.2);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==13){

            spEnergiesInterval->setMinE(1900.0);
            spEnergiesInterval->setMaxE(2227.82);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==14){

            spEnergiesInterval->setMinE(2227.82);
            spEnergiesInterval->setMaxE(2574.73);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==15){

            spEnergiesInterval->setMinE(2574.73);
            spEnergiesInterval->setMaxE(2972.15);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==16){

            spEnergiesInterval->setMinE(2972.15);
            spEnergiesInterval->setMaxE(3359.3);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==17){

            spEnergiesInterval->setMinE(3359.3);
            spEnergiesInterval->setMaxE(3772.04);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==18){

            spEnergiesInterval->setMinE(3772.04);
            spEnergiesInterval->setMaxE(4199.83);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==19){

            spEnergiesInterval->setMinE(4199.83);
            spEnergiesInterval->setMaxE(4648.26);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==20){

            spEnergiesInterval->setMinE(4648.26);
            spEnergiesInterval->setMaxE(5120.45);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==21){

            spEnergiesInterval->setMinE(5120.45);
            spEnergiesInterval->setMaxE(5598.72);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==22){

            spEnergiesInterval->setMinE(5598.72);
            spEnergiesInterval->setMaxE(6079.87);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==23){

            spEnergiesInterval->setMinE(6079.87);
            spEnergiesInterval->setMaxE(6724.18);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==24){

            spEnergiesInterval->setMinE(6724.18);
            spEnergiesInterval->setMaxE(7349);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==25){

            spEnergiesInterval->setMinE(7349);
            spEnergiesInterval->setMaxE(7961.68);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else if(z==26){

            spEnergiesInterval->setMinE(7961.68);
            spEnergiesInterval->setMaxE(9500);

            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));

        }
        else {
            spEnergiesInterval->setMinE(99999);
            spEnergiesInterval->setMaxE(99999);
            m_zVect.push_back(ZVectEntry2(z,spEnergiesInterval));
        }
    }//end of loop on atomic numbers


    //m_log << MSG::INFO << "GcrSelectTool END fillExpectedEdep()" << endreq ;  


}


//-----------------------------------------------------------------------------------------------------------------

void GcrSelectTool::verifyExpectedEdep (){

    //m_log << MSG::INFO << "GcrSelectTool BEGIN verifyExpectedEdep()" << endreq ;  

    double minE,maxE;



    for(ZVect::iterator zVecIter = m_zVect.begin(); zVecIter != m_zVect.end(); zVecIter++){//loop on atomic numbers
        ZVectEntry2 currentZVectEntry = *zVecIter;
        int currentZ = currentZVectEntry.getZ();
        SpEnergiesInterval* currentSpEnergiesInt = currentZVectEntry.getSpEnergiesInterval();
        minE = currentSpEnergiesInt->getMinE();
        maxE = currentSpEnergiesInt->getMaxE();
        m_log << MSG::INFO << "currentZ=" << currentZ <<  ", Emin=" << minE << ", Emax=" << maxE << endreq ;  

    }// end of loop on atomic numbers

    //m_log << MSG::INFO << "GcrSelectTool END verifyExpectedEdep()" << endreq ;  

}


//-----------------------------------------------------------------------------------------------------------------
int GcrSelectTool::inferZ(double correctedEnergy){

    m_log << MSG::INFO << "GcrSelectTool BEGIN inferZ(double)" << endreq ; 

    m_log << MSG::INFO << "correctedEnergy=" << correctedEnergy << endreq ;  

    double minE,maxE;
    int currentZ=-700;
    //m_log << MSG::INFO << "m_zVect.size()=" << m_zVect.size() << endreq ;  


    for(ZVect::iterator zVecIter = m_zVect.begin(); zVecIter != m_zVect.end(); zVecIter++){//loop on atomic numbers
        ZVectEntry2 currentZVectEntry = *zVecIter;
        currentZ = currentZVectEntry.getZ();
        //SpEnergiesVec* currentSpEnergiesVec = currentZVectEntry.getSpEnergiesVec();
        //LayerSpEnData firstLayerSpEnData = currentSpEnergiesVec->at(layerNum);
        //EPeak=firstLayerSpEnData.getEPeak();
        //sigma=firstLayerSpEnData.getSigma();

        SpEnergiesInterval* currentSpEnergiesInterval = currentZVectEntry.getSpEnergiesInterval();
        minE = currentSpEnergiesInterval->getMinE();
        maxE = currentSpEnergiesInterval->getMaxE();

        //m_log << MSG::INFO << "currentZ=" << currentZ << ", Emin=" << minE << ", Emax= " << maxE << endreq ;  


        //if( ((EPeak-n*sigma)<=correctedEnergy) && (correctedEnergy <= (EPeak+n*sigma)) ){
        if(  (minE<=correctedEnergy) && (correctedEnergy < maxE) ){
            m_log << MSG::INFO << "Zfound, Z=" << currentZ<< endreq ; 
            return currentZ;
        }       
    }// end of loop on atomic numbers



    m_log << MSG::INFO << "GcrSelectTool END inferZ(double)" << endreq ; 
    return -700;
    }


    //-----------------------------------------------------------------------------------------------------------------



    int GcrSelectTool::buildSelectedXtalsVec (){

        bool debugging = false;

        m_log << MSG::INFO << "GcrSelectTool BEGIN buildSelectedXtalsVec()" << endreq ;  

        m_gcrSelectedXtalsVec.clear();

        m_inferedZ = 0;

        int iitow,iilay,iicol;

        for(GcrTowersVec::iterator gcrTowersVecIter = m_gcrTowersVec.begin(); gcrTowersVecIter != m_gcrTowersVec.end(); gcrTowersVecIter++){
            GcrTower currentTow = *gcrTowersVecIter;
            GcrLayersVec& currentLayersVec = currentTow.getLayersVec();
            //m_log << MSG::INFO << "itow:"<< itow << " gcrLayersVec.size()= " << currentLayersVec->size() << endreq;

            bool firstGoodLayerFound = false;

            for(GcrLayersVec::iterator gcrLayersVecIter=
currentLayersVec.begin(); gcrLayersVecIter!= currentLayersVec.end(); gcrLayersVecIter++){

                GcrLayer currentLay = *gcrLayersVecIter;

                bool isGlobalLayerGood = ((m_layMatchTrackCrit[currentLay.getLayerNb()]==1)&& (m_layMatchMultCrit2[currentLay.getLayerNb()] == 1));

                if(debugging){
                    m_log << MSG::INFO << "itow,ilay:"<< currentTow.getTowerNb()<< "," << currentLay.getLayerNb() << endreq;
                    m_log << MSG::INFO << "is currentLay good?"<< currentLay.getIsGoodLayer() << endreq;
                    m_log << MSG::INFO << "is Global currentLay good?"<< isGlobalLayerGood<< endreq;
                }


                if(isGlobalLayerGood) // si la couche (globale) est bonne
                    //if(currentLay.getIsGoodLayer())
                    firstGoodLayerFound = true;
                else {  //si la couche (glosbale ou pour cette tour) n'est pas bonne       
                    if(!firstGoodLayerFound)  // et que la premiere couche bonne n'a pas ete trouvée, continue de chercher
                        continue;
                    else   // si la couche n'est pas bonne et que la premiere couche bonne a deja ete trouvée, on arrete tout
                        break;
                }


                GcrClustersVec* currentClustersVec = currentLay.getClustersVec();


                if(currentClustersVec != NULL){


                    for(GcrClustersVec::iterator currentClustersVecIter = currentClustersVec->begin(); currentClustersVecIter != currentClustersVec->end(); currentClustersVecIter++){

                        int cluLay = -1;
                        GcrCluster currentCluster = *currentClustersVecIter;
                        GcrHitsVec& currentHitsVec = currentCluster.getHitsVec();

                        if(currentHitsVec.size() > 0){

                            int currentHitNb = 0;

                            double rawEnergy;
                            double pathLength;
                            double correctedEnergy;

                            double closerFaceDist;
                            int crossedFaces;
                            Point entryPoint;
                            Point exitPoint;


                            for(GcrHitsVec::iterator currentHitsVecIter =
currentHitsVec.begin(); currentHitsVecIter != currentHitsVec.end(); currentHitsVecIter++){


                                GcrHit currentHit = *currentHitsVecIter;
                                //m_log << MSG::INFO << "currentHitNb="<< currentHitNb <<endreq;

                                if(currentHit.getIsGoodHit()){// for the moment, only Xtals that correspond exactly with theoretical track will be stored

                                    Event::CalXtalRecData* xtalData= currentHit.getXtalData();
                                    iitow=xtalData->getPackedId().getTower();
                                    iilay=xtalData->getPackedId().getLayer();
                                    iicol=xtalData->getPackedId().getColumn();

                                    cluLay=iilay;

                                    if(m_debugging)
                                        m_log << MSG::INFO << "grab GcrSelectedXtal["<< iitow << "," << iilay << "," << iicol << "]"<<endreq;

                                    Event::GcrXtal currentGcrXtal= m_gcrXtalVec.at(m_gcrXtalsMap[iitow][iilay][iicol]);
                                    rawEnergy = xtalData->getEnergy();// cette energie on l'a sur la table d'energies... cette table sera peut etre a enlever???
                                    pathLength    = currentGcrXtal.getPathLength();  /// ATTENTION: y a t'il la possibilite de pointer directement sur GcrXtal au lieu de sur CalXtalRecData dans GcrHit???, non(?) car on a besoin d'aller recuperer l'energie
                                    correctedEnergy = currentHit.getCorrectedEnergy();
                                    closerFaceDist=currentGcrXtal.getClosestFaceDist();
                                    crossedFaces  =currentGcrXtal.getCrossedFaces();
                                    entryPoint   = currentGcrXtal.getEntryPoint();
                                    exitPoint    = currentGcrXtal.getExitPoint();

                                    if(debugging)
                                        m_log << MSG::INFO << "corrEnergy= "<< correctedEnergy <<endreq;

                                    if (currentCluster.getTotalRawEnergy()>0){
                                        currentCluster.setTotalRawEnergy(currentCluster.getTotalRawEnergy()+rawEnergy);
                                        currentCluster.setTotalPathLength(currentCluster.getTotalPathLength()+pathLength);
                                    }
                                    else{
                                        currentCluster.setTotalRawEnergy(rawEnergy);
                                        currentCluster.setTotalPathLength(pathLength);

                                    }

                                    if(debugging){ 
                                        m_log << MSG::INFO << "Xtal["<< iitow << "," << iilay << "," << iicol <<"]:" 
                                            << "XtalpathLength= " << pathLength<< " XtalrawEnergy= " << rawEnergy
                                            << "currentCluster.getTotalRawEnergy()= " << currentCluster.getTotalRawEnergy()
                                            << "currentCluster.getTotalPathLength()= " << currentCluster.getTotalPathLength()
                                            << endreq;
                                    }


                                    if(debugging)
                                        m_log << MSG::INFO << "grab Xtal["<< iitow << "," << iilay << "," << iicol <<"]"<< endreq ; 

                                    m_gcrSelectedXtalsVec.push_back(Event::GcrSelectedXtal(xtalData->getPackedId(),rawEnergy,pathLength,correctedEnergy,-1,closerFaceDist,crossedFaces,entryPoint, exitPoint ));

                                }

                                //m_log << MSG::INFO << "XtalpackedId, iHit=" << currentHit.getXtalData()->getPackedId()<< ","<< currentHitNb << " currentHit corr energy= " << currentHit.getCorrectedEnergy() << endreq;
                                currentHitNb++;
                            }//end of iteration on GcrHitsVec

                        }// end of if(currentHitsVec)
                        else{
                            m_log << MSG::INFO << "no currentHitsVec" << endreq;
                        }

                        if(cluLay == 0){
                            if(debugging){
                                m_log << MSG::INFO << "Cluster TotalRawEnergy=" << currentCluster.getTotalRawEnergy()<< endreq;
                                m_log << MSG::INFO << "Cluster TotalPathLength=" << currentCluster.getTotalPathLength()<< endreq;
                            }
                            m_inferedZ = inferZ(m_xtalHeight*(currentCluster.getTotalRawEnergy() / currentCluster.getTotalPathLength()));

                        }




                    }/// end of iteration on clusters

                }///end of if(currentClustersVec)
                else{
                    //m_log << MSG::INFO << "no clustersVec found" << endreq;


                }



            }// end of loop on LayersVec


        }//end of loop on TowersVec


        m_log << MSG::INFO << "m_gcrSelectedXtalsVec.size()" << m_gcrSelectedXtalsVec.size()<< endreq;

        return m_gcrSelectedXtalsVec.size();


    }

    // ----------------------------------------------------------------------------

    /**
     * @author CL 06/02/2006
     * This method allows to store GcrSelectedXtals in TDS structure
     */
    StatusCode GcrSelectTool::storeGcrSelectedXtals () {
        StatusCode sc = StatusCode::SUCCESS;

        m_log << MSG::DEBUG << "BEGIN storeGcrSelectedXtals in GcrSelectTool" << endreq;

        m_gcrSelectedXtalsCol = SmartDataPtr<Event::GcrSelectedXtalsCol>(m_dataSvc,EventModel::CalRecon::GcrSelectedXtalsCol);

        // If no pointer then create it
        if (m_gcrSelectedXtalsCol == 0)
        {
            m_gcrSelectedXtalsCol = new Event::GcrSelectedXtalsCol();
            sc = m_dataSvc->registerObject(EventModel::CalRecon::GcrSelectedXtalsCol, m_gcrSelectedXtalsCol);
            if (sc.isFailure()) throw GaudiException("Failed to create GCR Selected Xtal Collection!", name(), sc);
        }

        int nbStoredSelGcrXtals=0;

        //m_log << MSG::INFO << "m_gcrSelectedXtalsVec.size()=" << m_gcrSelectedXtalsVec.size()  << endreq;



        for(Event::GcrSelectedXtalsVec::iterator GcrSelectedXtalIter=m_gcrSelectedXtalsVec.begin(); GcrSelectedXtalIter != m_gcrSelectedXtalsVec.end(); GcrSelectedXtalIter++)
        {
            Event::GcrSelectedXtal& gcrSelectedXtal = *GcrSelectedXtalIter;

            // Need to create a new GcrXtal which will be "owned" by the TDS
            Event::GcrSelectedXtal* newGcrSelectedXtal = new Event::GcrSelectedXtal();

            // Now copy to this new GcrXtal. This should properly copy all elements (including inherited vector)
            *newGcrSelectedXtal = gcrSelectedXtal;
            /**m_log << MSG::INFO << "Dans GcrSelectTool::storeGcrSelectedXtals, newGcrSelectedXtal->getInferedZ()= " << newGcrSelectedXtal->getInferedZ() 
              << "raw Energy=" << newGcrSelectedXtal->getRawEnergy()<< " corrEnergy= " << newGcrSelectedXtal->getCorrEnergy()
              << "entry=" << newGcrSelectedXtal->getEntryPoint() <<"exit=" << newGcrSelectedXtal->getExitPoint()
              << "crossedFaces="<< newGcrSelectedXtal->getCrossedFaces()  << "closerFaceDist="<< newGcrSelectedXtal->getCloserFaceDist()
              << endreq;*/

            // Store in collection (and reliquish ownership of the object)
            m_gcrSelectedXtalsCol->push_back(newGcrSelectedXtal); 
            nbStoredSelGcrXtals++;
        }

        m_log << MSG::INFO << "nbStoredSelGcrXtals=" << nbStoredSelGcrXtals  << endreq;

        return sc;


    }

    // ----------------------------------------------------------------------------

    /**
     * @author CL 06/02/2006
     * This method allows to store GCRSelectVals in TDS structure
     */
    StatusCode GcrSelectTool::storeGcrSelectVals () {
        StatusCode sc = StatusCode::SUCCESS;

        m_log << MSG::VERBOSE << "BEGIN storeGcrSelectVals in GcrSelectTool" << endreq;

        m_gcrSelectVals = SmartDataPtr<Event::GcrSelectVals>(m_dataSvc,EventModel::CalRecon::GcrSelectVals);

        // If no pointer then create it
        if (m_gcrSelectVals == 0)
        {
            m_gcrSelectVals = new Event::GcrSelectVals();
            sc = m_dataSvc->registerObject(EventModel::CalRecon::GcrSelectVals, m_gcrSelectVals);
            if (sc.isFailure()) throw GaudiException("Failed to create GCR SelectVals !", name(), sc);
        }


        if(m_debugging)
            m_log << MSG::INFO << "In storeGcrSelectVals, inferedZ =" << m_inferedZ << endreq;

        m_gcrSelectVals->setInferedZ(m_inferedZ);
        m_gcrSelectVals->setAcdZ(-999);
        m_gcrSelectVals->setInteractionParams(-999);

        //gcrOBFStatusWord comes from gcrReconVals - loaded from TDS in GcrSelectTool:chooseEth() method
        m_gcrSelectVals->setGcrOBFStatusWord(m_gcrOBFStatusWord);

        m_log << MSG::DEBUG << "In storeGcrSelectVals, gcrOBFStatusWord =" << m_gcrOBFStatusWord << endreq;

        m_log << MSG::VERBOSE << "END storeGcrSelectVals in GcrSelectTool" << endreq;

        return sc;


    }


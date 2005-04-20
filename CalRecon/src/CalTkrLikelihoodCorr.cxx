#include "CalTkrLikelihoodCorr.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(CalTkrLikelihoodCorr) ;

CalTkrLikelihoodCorr::CalTkrLikelihoodCorr( const std::string& type, const std::string& name, const IInterface* parent)
  :CalLikelihoodCorr(type,name,parent){
    // declare base interface for all consecutive concrete classes
    declareInterface<ICalEnergyCorr>(this);
    declareProperty("dataFile",
                    m_dataFile="$(CALRECONROOT)/xml/CalTkrLikelihood.data");
};

StatusCode CalTkrLikelihoodCorr::initialize(){
    
    // parent initialize must not be forgotten
    if (CalEnergyCorr::initialize().isFailure())
     { return StatusCode::FAILURE ; }

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	  log << MSG::INFO << "Initializing CalTkrLikelihoodCorr" <<endreq;

    typedef std::map<double*,std::string> PARAMAP;
    PARAMAP param;
    param[&m_pitchTOWER]= std::string("towerPitch");
    param[&m_halfwidthCAL]= std::string("CALModuleWidth");
    param[&m_heightCAL]= std::string("CALModuleHeight");

    // now try to find the GlastDevSvc service
    IGlastDetSvc* detSvc;
    sc= service("GlastDetSvc", detSvc);
    if (sc.isFailure() ) {
        log << MSG::ERROR << "  Unable to get GlastDetSvc " << endreq;
        return sc;
    }

    for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
      if(!detSvc->getNumericConstByName((*it).second, (*it).first)) {
        log<<MSG::ERROR<<" constant "<<(*it).second<<" not defined"<<endreq;
        return StatusCode::FAILURE;
      }
    }
    m_zOriginCAL= -47.395;
    m_halfwidthCAL*= .5;
    

    // Read in the parameters from the data file
    readPDFparameters( log, m_dataFile );
    return sc;
}


StatusCode CalTkrLikelihoodCorr::doEnergyCorr( Event::CalCluster * cluster )
//Purpose and method:
//
//   This function performs:
//   The main actions are:
//      - check wheter the event meets basic requirements (CUTS)
//      - calculate energy by TKR correction method using CalLikelihoodCorr 
// 
// TDS input: CalCluster (and TkrDigi?)
// TDS output: CalClusters
{
  StatusCode sc = StatusCode::SUCCESS;

  // if reconstructed tracker data doesn't exist or number of tracks is 0:
  if( !getKernel()->getTkrNVertices() ) {
    //cluster->setEnergyCalTkrLikelihood((double) CalLikelihoodCorr::cutNOTKRREC, 1.);
    return sc;
  }
      
  // CUTS
  // this checks whether a set of PDF parameters exist for this event's
  // energy, direction and point of impact.
  // if not, the setEnergyCorr value corresponds to the
  // defined CalTkrLikelihoodCorrFlags values
  
  // Energy:
  double pdfVariables[2]= {cluster->getEnergySum(), 0.};
  double pdfDataPoint[2]= {0., getKernel()->getSlope(cluster)};
  if( pdfVariables[0]<20. ) {
    //cluster->setEnergyCalTkrLikelihood((double) CalLikelihoodCorr::cutMAXCALENERGY, 1.);
    return sc;
  }
  if( pdfVariables[0]>1900. ) { 
    //cluster->setEnergyCalTkrLikelihood((double) CalLikelihoodCorr::cutMAXCALENERGY, 1.);
    return sc;
  }
  
  // direction: slope must be above \f$cos(32\circ)$\f
  if( getKernel()->getSlope(cluster)<.8480481 ){ 
    //cluster->setEnergyCalTkrLikelihood((double) CalLikelihoodCorr::cutSLOPE, 1.);
    return sc; 
  }

  // z-position:
  const Point &trackVertex= getKernel()->getTkrFrontVertex()->getPosition();
  int vertex=  int((trackVertex[2]-108.)*3.2e-2);
  if( vertex<0 || vertex>15 ) { 
    //cluster->setEnergyCalTkrLikelihood((double) CalLikelihoodCorr::cutVERTEX, 1.);
    return sc; 
  }

  // x,y-position: 1. find the track position on top of the CAL
  //               2. calculate its distance to parallel sides of the CAL,
  //    integrated along the track, normalized to 1.
  //               3. use the geometric mean of values on X and Y axis as a
  //    basis for geometric cut.
  const Vector &trackDir= getKernel()->getTkrFrontVertex()->getDirection();
  double geometricCut= 1.;
  for( int ax= 0; ax<2; ++ax ){
    double slope= -trackDir[ax]/trackDir[2];
    double posCAL= trackVertex[ax];
    posCAL+= (trackVertex[2]-m_zOriginCAL)*slope;
    geometricCut*= integratedDistance2TowerSide(posCAL, slope);
  }
  geometricCut= sqrt(geometricCut);

  // correcting for event slope
  double thetaNorm= .707106781186547462*sqrt(1./(trackDir[2]*trackDir[2])-1.);
  thetaNorm*= m_heightCAL/(2*m_halfwidthCAL); 
  geometricCut/= (1-thetaNorm);
  
  if( geometricCut<.01 ){
    //cluster->setEnergyCalTkrLikelihood((double) cutPOSITION, 1.);
    return sc;
  }
  int geometricCutindex= (geometricCut>.35)+(geometricCut>.75);

  // CALCULATE
  //    - get number of hits in TKR
  Event::TkrDigiCol *tkrDigiData = SmartDataPtr<Event::TkrDigiCol>
                          (getKernel()->getEventSvc(), EventModel::Digi::TkrDigiCol); 
  Event::TkrDigiCol::iterator it;
  for ( it= tkrDigiData->begin(); it!=tkrDigiData->end(); ++it )
    pdfVariables[0]+= (double)(*it)->getNumHits();

  double recEnergy= 0., recEnergyWidth= 1.;

  sc= calculateEvent( vertex+geometricCutindex*15,
                      50., 2000., pdfVariables, pdfDataPoint,
                      recEnergy, recEnergyWidth );

  //cluster->setEnergyCalTkrLikelihood(recEnergy, recEnergyWidth);
  return sc;
}

double CalTkrLikelihoodCorr::integratedDistance2TowerSide( double x, 
                                                  double slope ) const {
  if( slope==0. ) return 1.;
  if( slope>0. ){ x*= -1.; slope*= -1; }
  bool leaving= false;
  if( x<0. ) x+= m_pitchTOWER;
  if( x<0. ) x+= m_pitchTOWER;
  if( x>m_pitchTOWER ){ x-= m_pitchTOWER; leaving= true; }
  x-= m_pitchTOWER*.5;
  if( fabs(x)>m_pitchTOWER*.5 ) return 0.;
  if( x>m_halfwidthCAL ){
    if( leaving ) return 0.;
    x-= m_pitchTOWER;
  }
  double height= m_heightCAL;
  if( x<-m_halfwidthCAL ) {
    if( (height-= (m_halfwidthCAL+x)/slope)<=0. ) return 0.;
    x= -m_halfwidthCAL;
  } 
  double result= 0.;
  double norm= height*m_halfwidthCAL;

  if( x<0 ){
    double z= -x/slope;
    if( -z>=height ) {
      result= (m_halfwidthCAL+x-.5*height*slope)*height/norm;
      return result;
    }
    result= -(2.*m_halfwidthCAL+x)*z*.5;
    height+= z;
    x= 0.;
  } 
  double z= (m_halfwidthCAL-x)/slope;
  if( -z>=height ) result+= (m_halfwidthCAL-x+height*slope*.5)*height;
  else  result-= (m_halfwidthCAL-x)*z*.5;
  result/= norm;
  return result;
}

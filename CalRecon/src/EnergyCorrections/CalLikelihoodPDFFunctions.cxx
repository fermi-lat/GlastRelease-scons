
#include "GaudiKernel/MsgStream.h"
#include "CalLikelihoodPDFFunctions.h"
#include "CalLikelihoodManagerTool.h"
#include "Event/Digi/TkrDigi.h"

/******************************************************************************/
/***********************  PDFFunction *****************************************/
/******************************************************************************/
void PDFFunction::clear(void){
  if( m_NintPar ) delete m_IntPar;
  if( m_NfcnPar ) delete m_FcnPar;
  if( m_NfcnArg ) delete m_FcnArg;
}
void PDFFunction::setNfunctionParameters(int n){
  m_NfcnPar= n;
  if( m_FcnPar ) delete []m_FcnPar;
  m_FcnPar= n?new double[n]:0;
}
void PDFFunction::setNfunctionArgumentents(int n){
  m_NfcnArg= n;
  if( m_FcnArg ) delete []m_FcnArg;
  m_FcnArg= n?new double[n]:0;
}
void PDFFunction::setNinterpolationParameters(int n){
  m_NintPar= n;
  if( m_IntPar ) delete []m_IntPar;
  m_IntPar= n?new double[n]:0;
}

bool PDFFunction::eval(PDFParameters *table, double ans[1], bool init)
{ return table->interpolation(m_IntPar, m_FcnPar, init) || value(ans); }

bool PDFFunction::getMPV(PDFParameters* table, double eMin, double eMax,
                         double mpv[2]){
  // arguments are:
  //  - eMin, eMax: boundaries for search

  // see CalRecon CalLikelihoodTool for information 
  // next values are the range boundaries
  trialEnergy()= eMin;
  if( table->initialise(m_IntPar) ) return true;
  
  // next variables are for the estimation of the quality of the 
  // reconstruction: m_widthEnergyCorr
  double &eRecon= mpv[0];
  double &probRecon= mpv[1];
  eRecon= eMin;
  probRecon= -1.e40;
  for( double bW= delta(eMin, eMax)/5; bW>.1; bW= delta(eMin, eMax) )
  {
    int errCalls= 0;
    for( trialEnergy()= eMin; trialEnergy()<eMax; trialEnergy()+= bW )
    {
      // get log normal parameters
      double probTrial;
      if( eval(table, &probTrial, false) ) 
      { 
        ++errCalls; 
        continue; 
      }

      // calculate pdf value with these parameters
      if( probTrial>probRecon ) 
      { 
        probRecon= probTrial; 
        eRecon= trialEnergy();
      }
    }
    if( errCalls==m_Nstep ) return true;

    eMin= eRecon-bW;
    eMax= eRecon+bW;
    if( eMin<table->minTrialEnergy() ) eMin= table->minTrialEnergy();
    if( eMax>table->maxTrialEnergy() ) eMax= table->maxTrialEnergy();
  }

  return fabs(eRecon-table->minTrialEnergy())<1. 
      || fabs(eRecon-table->maxTrialEnergy())<1.;
}

bool PDFFunction::getFWHM(PDFParameters *table, const double mpv[2], 
                                           double fwhmLimits[2] ){
  // QUALITY
  // now looking for FWHM
  // fwhmLimits is an outer limit for an energy such that 
  // eval(x)<.5*maxProb
  if( table->initialise(m_IntPar) ) return true;
  fwhmLimits[0]= table->minTrialEnergy();
  fwhmLimits[1]= table->maxTrialEnergy();
  double fwhmValues[2]= {0., 0.};

  for( int iFWHM= 0; iFWHM<2; ++iFWHM )
  {
    // find starting point
    trialEnergy()= mpv[0];
    do {
      trialEnergy()*= iFWHM?1.1:.9;
      if( trialEnergy()>fwhmLimits[iFWHM] ) {
        if( eval(table, fwhmValues+iFWHM, false) ) fwhmValues[iFWHM]= 0.;
        fwhmValues[iFWHM]/= mpv[1];
      }
    } while ( trialEnergy()>fwhmLimits[iFWHM] && fwhmValues[iFWHM]>.5 );

    if( fwhmValues[iFWHM]>.5 ) { 
      // in case of no starting point, move first arg around until 1 is found
      double saving1stArg= firstArg();
      double calE[2]= {firstArg()*(iFWHM?.2:1.), iFWHM?firstArg():mpv[0]};
      double newMPV[2]= {mpv[0], mpv[1]};
      if( calE[1]<table->getGrid()->getBinCenter(0, 3)*.2 ) 
        calE[1]= table->getGrid()->getBinCenter(0, 3);
      while( fabs(fwhmValues[iFWHM]-.5)>5e-3 && fabs(calE[0]-calE[1])>.1 ){
        firstArg()= (calE[0]+calE[1])*.5;
        table->initialise(m_IntPar);
        if( getMPV(table, firstArg(), firstArg()*5, newMPV) ) calE[0]= calE[1];
        else {
          trialEnergy()= fwhmLimits[iFWHM];
          if( eval(table, fwhmValues+iFWHM, false) 
              || (fwhmValues[iFWHM]/= newMPV[1])>1. ) calE[0]= calE[1];
          else calE[iFWHM ^ (fwhmValues[iFWHM]<.5)]= firstArg();
        }
      }
      firstArg()= saving1stArg;
      if( fabs(fwhmValues[iFWHM]-.5)>.01 ) fwhmLimits[iFWHM]= -1.;
      else fwhmLimits[iFWHM]= fabs(1.-fwhmLimits[iFWHM]/newMPV[0]);
    } else {
      double lim[2]= {iFWHM?mpv[0]:trialEnergy(), iFWHM?trialEnergy():mpv[0]};
      int errCalls= 0;
      for( double bW= delta(lim[0], lim[1]); bW>.1; bW= delta(lim[0], lim[1]) )
      {
        errCalls= 0;
        double maxE= lim[1];
        for( trialEnergy()= lim[0]; trialEnergy()<maxE; trialEnergy()+= bW )
        {
          double probTrial;
          if( eval(table, &probTrial, false) )
          { 
            ++errCalls; 
            continue;
          } 
          if( (probTrial<mpv[1]*.5) && (iFWHM^(trialEnergy()>lim[iFWHM])) )
            lim[iFWHM]= trialEnergy();
        }
        if( errCalls==m_Nstep )  {
          fwhmLimits[iFWHM]= -1.;
          break;
        }
        lim[!iFWHM]= lim[iFWHM]+(1-2*iFWHM)*bW;
        lim[iFWHM]-= (1-2*iFWHM)*bW;
      }
      fwhmLimits[iFWHM]= fabs(1-lim[iFWHM]/mpv[0]);
    }
  }

  switch ( (fwhmLimits[0]==-1.) + 2*(fwhmLimits[1]==-1.) ) {
    case 3:
      return true;
    case 2:
      fwhmLimits[1]= fwhmLimits[0];
      return false; 
    case 1:
      fwhmLimits[0]= fwhmLimits[1];
      return false;
  } 
  return false;
}
/******************************************************************************/
/***********************  PDFLikelihood ***************************************/
/******************************************************************************/
void PDFLikelihood::setEvt(const Event::CalCluster*cluster,
                           const Event::TkrVertex*vertex) {
  tkr1ZDir()= fabs(vertex->getDirection()[2]);
  calEnergyRaw()= cluster->getCalParams().getEnergy();
  calELayer7()= cluster->back().getEnergy();

  const Event::TkrDigiCol *digi= m_manager->getTkrDigiCol();
  int nHits= 0;
  for(Event::TkrDigiCol::const_iterator hit= digi->begin(); hit!=digi->end();
       ++hit)
    nHits+= (*hit)->getNumHits();
  tkrSumHits()= nHits;
}
bool PDFLikelihood::value(double result[1]) const{
  // calculate LogNormal(x):
  // \f[ LogNormal(x, parameters(recEnergy)) \f]
  // if( par[3]==0. ) return 0.; that shouldn't happen
  // \f[ LogNormal(x)= N\,\exp(-\frac{1}{2}(\frac{\ln(1+\tau (x-\mu) \frac{\sinh({\tau}\sqrt{\ln\,4})}{ 2.36\beta{\tau}\,\sqrt{\ln\,4}})}\tau)^2+\tau^2) \f]

  //  parameter (TKR Hits or CAL Last Layer Energy)
  if( fabs(sigma())<1.e-10 ) return true;

  errno = 0 ;
  double shTau= tau()*1.17741002251547466;  //sqrt(log(4))
  result[0]= calEnergyRaw()+calELayer7()*calELayer7Alpha();
  result[0]+= tkrSumHits()*tkrSumHitsAlpha()-mpv();
  result[0]/= sigma();

  if( fabs(shTau)>1.e-10 )
  {
    shTau= sinh(shTau)/(shTau);
    result[0]= log(1+tau()*result[0]*shTau)/tau(); 
  } 
  else shTau= 1.;
  
  result[0]= exp(-0.5*(result[0]*result[0]+tau()*tau()));

  // sqrt(2*pi)
  result[0]*= occupancy()*norm()*shTau/(2.50662827463100024*sigma());

  if( errno ) result[0]= 0.;
  return false;
}
/******************************************************************************/
/***********************  PDFLowEnergyCuts ************************************/
/******************************************************************************/
PDFLowEnergyCuts::PDFLowEnergyCuts(const CalLikelihoodManagerTool*manager,
                                   MsgStream &log )
    : PDFFunction(2, 2, 2){
  typedef std::map<double*,std::string> PARAMAP;
  PARAMAP param;
  param[&m_towerPitch]  = std::string("towerPitch");
  param[&m_ratioCDEHeighTowerPitch]= std::string("cellVertPitch");
  manager->getParameters(param, log);
  m_calZorigin    = -47.395;
  m_ratioCDEHeighTowerPitch /= m_towerPitch;
}

bool PDFLowEnergyCuts::value(double result[1]) const{
  result[0]= -1.;
  if( calZEcntr()<getFunctionArguments()[cZECNTR_MIN] ) return true;
  if( calZEcntr()>getFunctionArguments()[cZECNTR_MAX] ) return true;
  if( geometricCut()<.17 ) return true;
  if( geometricCut()>.60 ) result[0]= 3.;
  else if( geometricCut()>.40 ) result[0]= 2.;
  else if( geometricCut()>.25 ) result[0]= 1.;
  else result[0]= 0.;
  return false;
}
void PDFLowEnergyCuts::setEvt(const Event::CalCluster*cluster,
                              const Event::TkrVertex*vertex) {
  tkr1ZDir()= fabs(vertex->getDirection()[2]);
  calZEcntr()= cluster->getCalParams().getCentroid().z();
  geometricCut()= geometricCut(cluster, vertex);
}
double PDFLowEnergyCuts::geometricCut(const Event::CalCluster*cluster,
                                      const Event::TkrVertex*vertex) const
// finds the vertex position and a value equal to the integration of the
// distance to the closest crack along the trajectory, weighted by the energy
// in the layer.
{
  Vector pX= vertex->getDirection(); 
  Point x= vertex->getPosition(); 

  if( fabs(pX[2])<1e-10 ) return 0;
  double geometricCut= 0.;
  double weight= 0.;

  double slopes[2]= { pX[0]/pX[2], 
                      pX[1]/pX[2] };
  // next values are for normalisation purposes: 
  // they correct differences to a normal incident particle
  double ellCorr[2]= {sqrt(1+slopes[0]*slopes[0]),
                      sqrt(1+slopes[1]*slopes[1])};
  double xT[2]= {x[0], 
                 x[1]};
  for( int ax=0; ax<2; ++ax )
  {
    // m_calZorigin: origin fixed at -47.395.
    // XY position is extrapolated to there
    xT[ax]-= slopes[ax]*(x[2]-m_calZorigin);
    // m_towerPitch: tower width
    // xT is in units of tower width
    xT[ax]/= m_towerPitch;
  
    // this is one step along the trajectory on the ax axis
    // There are 10 steps in all: thus the .1 * (CDE height)/(TOWER width)
    slopes[ax]*= .1*m_ratioCDEHeighTowerPitch;
  }
  
  for(Event::CalCluster::const_iterator layer= cluster->begin();
      layer!=cluster->end(); ++layer ){
    weight+= (*layer).getEnergy();
    if( (fabs(xT[0])>2.) || (fabs(xT[1])>2.) ) continue;
    double val= 0.;
    for( int ii=0; ii<10; ++ii )
    {
      double towerX[2]= { xT[0], xT[1] };
      for( int ax=0; ax<2; ++ax )
      {
        // find position relative to tower center
        // the tower being whichever one xT is now at
        if( towerX[ax]<-1. ) towerX[ax]= fabs(towerX[ax]+1.5);
        else if( towerX[ax]<0. ) towerX[ax]= fabs(towerX[ax]+.5);
        else if( towerX[ax]<1. ) towerX[ax]= fabs(towerX[ax]-.5);
        else towerX[ax]= fabs(towerX[ax]-1.5);

        // normalisation for the trajectory slant
        towerX[ax]= (.5-(towerX[ax]<.5?towerX[ax]:.5))/ellCorr[ax];

        // moving along the trajectory
        xT[ax]-= slopes[ax];
      }

      // adding minimum distance to the tower edge

      val+= towerX[towerX[0]>towerX[1]];
    }
    geometricCut+= val*(*layer).getEnergy();
  }
  geometricCut*= .2*sqrt(.5/(pX[2]*pX[2])+.5)/weight;
  return geometricCut;
}

/******************************************************************************/
/***********************  PDFHighEnergyCuts ***********************************/
/******************************************************************************/
PDFHighEnergyCuts::PDFHighEnergyCuts(const CalLikelihoodManagerTool*manager,
                                     MsgStream &log )
    : PDFFunction(2, 7, 3){
  std::map<double*,std::string> param;
  param[&m_towerPitch]  = std::string("towerPitch");

  manager->getParameters(param, log);
}
bool PDFHighEnergyCuts::value(double result[1]) const{
  result[0]-= 1.;
  
  const double *cuts= getFunctionParameters();
  double bottom_out= 3*(-220)-2*cuts[cZECNTR_BOTTOM];
  double core_out= 1.5*cuts[cZECNTR_BOTTOM]-.5*cuts[cZECNTR_CORE];
  double top_out   = 3*cuts[cZECNTR_TOP]    - 2*cuts[cZECNTR_CORE];
  if( calZEcntr()>top_out ) return true;
  if( calZEcntr()<bottom_out ) return true;
  if( calTwrEdgeCntr()>cuts[cTOWER_CENTER] && calELayer7()>cuts[cE7] ) {
    if( calZEcntr()<core_out ) result[0]= 1.;
    else result[0]= 3.; 
    return false;
  }

  if( calELayer7()>cuts[cE7] ){
    if( calTwrEdgeCntr()>cuts[cTOWER_BORDER] )
      result[0]= 1.+(calZEcntr()>cuts[cZECNTR_BOTTOM]);
    
    else if( calTwrEdgeCntr()>cuts[cTOWER_EDGE] )
      result[0]= 1.+(calZEcntr()>cuts[cZECNTR_CORE]);

    else result[0]= (calZEcntr()<cuts[cZECNTR_TOP]);
    return false;
  }
  if( calZEcntr()>cuts[cZECNTR_BOTTOM]) result[0]= 0.;
  else return true;
  return false;
}
void PDFHighEnergyCuts::setEvt(const Event::CalCluster*cluster,
                              const Event::TkrVertex*vertex) {
  tkr1ZDir()= fabs(vertex->getDirection()[2]);
  calZEcntr()= cluster->getCalParams().getCentroid().z();
  calTwrEdgeCntr()= calTwrEdgeCntr(cluster);
  calELayer7()= cluster->back().getEnergy();
}
double PDFHighEnergyCuts::calTwrEdgeCntr(const Event::CalCluster *cluster)const{
  double x= cluster->getCalParams().getCentroid().x();
  double y= cluster->getCalParams().getCentroid().x();

  if(x<0) x+= m_towerPitch;
  if(x>m_towerPitch) x-= m_towerPitch;
  x= m_towerPitch-fabs(x-m_towerPitch*.5);
  
  if(y<0) x+= m_towerPitch;
  if(y>m_towerPitch) y-= m_towerPitch;
  y= m_towerPitch-fabs(y-m_towerPitch*.5);

  return x<y?x:y;
}


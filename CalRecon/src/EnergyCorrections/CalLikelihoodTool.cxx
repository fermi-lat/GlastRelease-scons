#include "CalLikelihoodTool.h"
#include "GaudiKernel/MsgStream.h" 
#include "GaudiKernel/GaudiException.h" 
#include "facilities/Util.h"
#include <fstream>

StatusCode CalLikelihoodTool::initialize()
{
    // This function does following initialization actions:
    //    - extracts geometry constants from xml file using GlastDetSvc
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
      log << MSG::INFO << "Initializing CalLikelihoodTool" <<endreq;

    typedef std::map<double*,std::string> PARAMAP;
    PARAMAP param;
    param[&m_towerPitch]  = std::string("towerPitch");
    param[&m_ratioCDEHeighTowerPitch]= std::string("cellVertPitch");

    //Locate and store a pointer to the data service which allows access to the TDS
    if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    if ((sc = service("GlastDetSvc", m_detSvc, true)).isFailure())
    { 
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }

    for(PARAMAP::iterator it=param.begin(); it!=param.end();it++)
    {
        if(!m_detSvc->getNumericConstByName((*it).second, (*it).first)) 
        {
            log<<MSG::ERROR<<" constant "<<(*it).second<<" not defined"<<endreq;
            throw GaudiException("Bad constant name ", name(), StatusCode::FAILURE);
        }
    }

    m_calZorigin    = -47.395;
    m_ratioCDEHeighTowerPitch /= m_towerPitch;
    
    // This function does following initialization actions:
    //    - extracts PDF parameters from txt file
	  log<< MSG::DEBUG <<"Initializing CalLikelihoodTool with file "
                     <<m_dataFile<<endreq;

    // Read in the parameters from the data file
    facilities::Util::expandEnvVar(&m_dataFile);
    std::ifstream dataFile(m_dataFile.data());
    if( !dataFile ){
      log<<MSG::ERROR<<" Unable to open data file: "<<m_dataFile<<endreq;
      return StatusCode::FAILURE;
    }
    
    int nPar;
    std::string line;
    getline(dataFile, line);

    m_Axes= 0;
    m_Data= 0;
    int pdf= 0;
    for( int lineNumber= 0; !dataFile.eof(); ++lineNumber )
    {
      getline(dataFile, line);
      if( line.size()>6 && (line.substr(0, 7)=="COMMENT") )
        continue; 

      else if( line.size()>6 && (line.substr(0, 6)=="#PDFs:") )
      {
        sscanf(line.data(), "#PDFs: %d #Parameters: %d", &m_Npdf, &nPar);
        log<<MSG::DEBUG<<"File must contain "<<m_Npdf<<", for "<<nPar
                       <<" parameters"<<endreq;
        m_Data= new PDF_Data*[m_Npdf];
        for( int pdfs= 0; pdfs<m_Npdf; ++pdfs ) 
          m_Data[pdfs]= 0;
        if( nPar>0 && m_Npdf>0 ) continue;
      } 

      else if( line=="AXES:" )
      {
        if( !m_Axes )
        {
          log<<MSG::DEBUG<<"Extracting Axes"<<endreq;
          if( (m_Axes= PDF_Axes::read(dataFile)) )
            continue;
          log<<MSG::ERROR<<"Axes did not build correctly"<<endreq;
        }
      }
      
      else if( line.size()>4 && line.substr(0,4)=="PDF:" )
      {
        if( !m_Data )
        {
          log<<MSG::ERROR<<"Line \"#PDFs: <number of PDFs> "
                           "Parameters: <number of parameters for PDF function>"
                           "\" must come  before 1st PDF";
        } 
        else if( pdf<m_Npdf )
        {
          log<<MSG::DEBUG<<"Extracting PDF Parameters for PDF #"<<pdf<<endreq;
          if( (m_Data[pdf++]= PDF_Data::read(dataFile, nPar, m_Axes)) )
              continue;
          log<<MSG::ERROR<<"PDF #"<<pdf-1<<" did not build correctly"<<endreq;
        }
      } 
      else if ( (line=="") || (line=="\n") 
                           || strspn(line.data(), " ")==strlen(line.data()) )
        continue;
            
      log<<MSG::ERROR<<"Incorrect data in file: "<<m_dataFile<<endreq
                     <<"Stopping on line :\""<<line.data()<<endreq;
      dataFile.close();
      return StatusCode::FAILURE;
  }

  dataFile.close();
  if( pdf<m_Npdf || !m_Axes )
  {
    log<<MSG::ERROR<<"Missing data in file: "<<m_dataFile<<endreq;
    return StatusCode::FAILURE;
  }
  return sc;
}

StatusCode CalLikelihoodTool::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    if( m_Data )
    {
        for( int plane= 0; plane<m_Npdf; ++plane )
        {
            if( m_Data[plane] ){ delete m_Data[plane]; m_Data[plane]= 0; }
        }
        delete []m_Data; 
    }
  
    if( m_Axes ) delete m_Axes;
  
    return sc;
}

Event::CalCorToolResult* CalLikelihoodTool::calculateEvent(
                                            const Event::CalCluster* cluster,
                                            MsgStream &log )
//Purpose and method:
//   This function performs:
//        Evaluation of the energy for a set of parameters (m_eventPar),
//      specifying the form of a function (evaluatePDF) which is evaluated at a
//      given point (pdfVariables, reconstruction energy). 
//        The function parameters are evaluated using linear in interpolation
//      of data points using m_eventPDF. 
//        The energy is estimated as the PDF's MPV for that given set of data
//      and variables. The quality of the reconstruction is estimated as the
//      FWHM of the PDF.
//        Both values are estimated recusively, calculating 10 points and
//      choosing a range around that turn's best estimate.
//      
//   The main actions are:
//      - estimate particle energy by TKR correction method (CALCULATE)
//      - estimate FWHM of the PDF for set a of parameters: 
//          m_eventPar && reconstructed energy            (QUALITY)
// 
// TDS input:
// TDS output:
{
  // CalCulate
  // Evaluation of energy using correlation method
  // Coefficients fitted using GlastSim.
  // looking for minimum by:
  //    -calculating 10 values on a range [x0, x1]
  //    -choosing a new range [mpv-binWidth, mpv+binWidth]
  //  and so on and so forth
  // every call, the range is decreased to less than range*.3
  // => max number of calls<log(MAXENERGY)/log(3) : calls to pdf fcn<100

  // pdf function is:
  // both the LogNormal parameters, mpv, width,..., and $\alpha$ applied
  // to nTkrHits depend on event vertex, direction and recEnergy

  // next values are the range boundaries
  log << MSG::DEBUG << "Starting Calculations" <<endreq;
  double mpv[2], fwhm[2];
  if( getMPV(mpv) || getFWHM(mpv, fwhm) ) return 0;
  log << MSG::DEBUG << "End of Calculations:"<< endreq 
      << MSG::VERBOSE << "Found Energy:" << mpv[0] 
                      << "MeV,\tWidth:" << .5*(fwhm[0]+fwhm[1])
                      << "\%" << endreq;


  // Create a CalCorToolResult object to hold the information
  Event::CalCorToolResult *corResult = new Event::CalCorToolResult();

  corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
  corResult->setCorrectionName(type());
  corResult->setChiSquare(1.);

  Event::CalParams params= cluster->getCalParams();
  params.setEnergy(mpv[0]);
  params.setEnergyErr(.5*(fwhm[0]+fwhm[1]));
  corResult->setParams(params);
  return corResult;
}

bool CalLikelihoodTool::getMPV(double mpv[2]){
  // next values are the range boundaries
  double lim[2]= { calEnergy(), calEnergy()*5.};
  
  // when this happens, must not set the bioudary to high or too low
  // too high and the steps will be too big
  if( lim[1] < m_Axes->getBinCenter(0, 4) )
    lim[1]= m_Axes->getBinCenter(0, 4);
  if( lim[0] < minTrialEnergy() ) lim[0] =  minTrialEnergy();
  if( lim[1] > maxTrialEnergy() ) lim[1] =  maxTrialEnergy();

  // next variables are for the estimation of the quality of the 
  // reconstruction: m_widthEnergyCorr
  double recEnergy= lim[0];
  double maxProb= -1e40;

  trialEnergy()= lim[0];
  m_Axes->init(m_eventPar);
  const int nSteps= 15;
  for( double bW= (lim[1]-lim[0])/nSteps; bW>.1; bW= (lim[1]-lim[0])/nSteps )
  {
    int errCalls= 0;
    for( trialEnergy()= lim[0]; trialEnergy()<lim[1]; trialEnergy()+= bW )
    {
      // get log normal parameters
      double trialProb;
      if( evaluatePDF(trialProb) ) 
      { 
        ++errCalls; 
        continue; 
      }

      // calculate pdf value with these parameters
      if( trialProb>maxProb ) 
      { 
        maxProb= trialProb; 
        recEnergy= trialEnergy();
      }
    }
    if( errCalls==nSteps ) return true;

    lim[0]= recEnergy-bW;
    lim[1]= recEnergy+bW;
    if( lim[0]<minTrialEnergy() ) lim[0]= minTrialEnergy();
    if( lim[1]>maxTrialEnergy() ) lim[1]= maxTrialEnergy();
  }
  mpv[1]= maxProb;
  mpv[0]= recEnergy;
  return fabs(mpv[0]-minTrialEnergy())<1. || fabs(mpv[0]-maxTrialEnergy())<1.; 
}
bool CalLikelihoodTool::getFWHM(const double mpv[2], double fwhmLimits[2]){
  // QUALITY
  // now looking for FWHM
  // fwhmLimits is an outer limit for an energy such that 
  // evaluatePDF(x)<.5*maxProb
  fwhmLimits[0]= minTrialEnergy()+.01;
  fwhmLimits[1]= maxTrialEnergy()-.01;
  const int nSteps= 15;
  for( int iFWHM= 0; iFWHM<2; ++iFWHM )
  {
    double delta= 0;
    trialEnergy()= fwhmLimits[iFWHM];
    if( evaluatePDF(delta) ) delta= 0.;
    delta/= mpv[1];
    if( delta<.5 ) 
    {
      double lim[2]= {iFWHM?mpv[0]:minTrialEnergy(),
                        iFWHM?maxTrialEnergy():mpv[0]};
      int errCalls= 0;
      for( double bW= (lim[1]-lim[0])/nSteps; bW>.1;bW= (lim[1]-lim[0])/nSteps )
      {
        errCalls= 0;
        double maxE= lim[1];
        for( trialEnergy()= lim[0]; trialEnergy()<maxE; trialEnergy()+= bW )
        {
          double trialProb;
          if( evaluatePDF(trialProb) )
          { 
            ++errCalls; 
            continue;
          } 
          // for lower(upper) FWHM x axis value, move the lim up(down) when
          // at a y axis value  below .5*maximum
          if( (trialProb<mpv[1]*.5) && (iFWHM^(trialEnergy()>lim[iFWHM])) )
            lim[iFWHM]= trialEnergy();
        }
        if( errCalls==nSteps ) return true;
        lim[!iFWHM]= lim[iFWHM]+(1-2*iFWHM)*bW;
        lim[iFWHM]-= (1-2*iFWHM)*bW;
      }
      fwhmLimits[iFWHM]= fabs(1-lim[iFWHM]/mpv[0]);
    } 
    else 
    {
      double savingCalE= calEnergy();
      double calE[2]= {calEnergy()*(iFWHM?.3:1.), calEnergy()*(iFWHM?1.:5.)};
      double newMPV[2];
      // when the lower(upper) FWHM x axis value is outside the PDF pahse
      // space, estimate it as close as possible from that point:
      // we move calEnergy() energy around to find it
      if( calE[1]<.1 ) calE[1]= m_Axes->getBinCenter(0, 4);
      while( fabs(delta-.5)>1e-2 && fabs(calE[0]-calE[1])>.1 )
      {
        calEnergy()= (calE[0]+calE[1])*.5;
        m_Axes->init(m_eventPar);
        if( getMPV(newMPV) ) return true; 
        trialEnergy()= fwhmLimits[iFWHM];
        if( evaluatePDF(delta) || (delta/= newMPV[1])>1. ) return true;
        calE[iFWHM ^ (delta<.5)]= calEnergy();
      }
      if( fabs(delta-.5)>.01 ) return true; 
      fwhmLimits[iFWHM]= fabs(1-fwhmLimits[iFWHM]/newMPV[0]);
      calEnergy()= savingCalE;
    }
  }
  return fwhmLimits[0]+fwhmLimits[1]>2.;
}

double CalLikelihoodTool::findGeometricCut( const Point &x,
                                            const Vector &pX, 
                                            const Event::CalCluster* cluster)
// finds the vertex position and a value equal to the integration of the
// distance to the closest crack along the trajectory, weighted by the energy
// in the layer.
// it is then translated in a cut value equal to the index of the PDF in the
// PDF table.
// this method will also be called by the tool's calibration alg.
{
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
  
  for ( Event::CalCluster::const_iterator layer= cluster->begin(); 
        layer!=cluster->end(); ++layer )
  {
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

bool CalLikelihoodTool::evaluatePDF( double &result ) const{
  // calculate LogNormal(x):
  // \f[ LogNormal(x, parameters(recEnergy)) \f]
  // if( par[3]==0. ) return 0.; that shouldn't happen
  // \f[ LogNormal(x)= N\,\exp(-\frac{1}{2}(\frac{\ln(1+\tau (x-\mu) \frac{\sinh({\tau}\sqrt{\ln\,4})}{ 2.36\beta{\tau}\,\sqrt{\ln\,4}})}\tau)^2+\tau^2) \f]

  double values[5]= {0., 0., 0., 0., 0. };
  //  parameter (TKR Hits or CAL Last Layer Energy)
  if( m_eventPDF->evaluateParameters(m_eventPar, values+1-hasCorrelation()) )
    return true;
  if( fabs(values[2])<1.e-10 ) return true;

  errno = 0 ;
  double shTau= values[4]*1.17741002251547466;  //sqrt(log(4))
  result= (m_eventPar[2]+m_eventPar[3]*values[0]-values[2])/values[3];

  if( fabs(values[4])>1.e-10 )
  {
    shTau= sinh(shTau)/(shTau);
    result= log(1+values[4]*result*shTau)/values[4]; 
  } 
  else shTau= 1.;
  
  result= exp(-0.5*(result*result+values[4]*values[4]));
  result*= values[1]*shTau/(2.50662827463100024*values[3]); //sqrt(2*pi)

  if( errno ) result= 0.;
  return false;
}


PDF_Axes *PDF_Axes::read(std::ifstream &dataFile)
{
  int nAxes, *sizes;
  dataFile>>nAxes;
  sizes= new int[nAxes];

  // extract axis sizes
  int ax= 0;
  int nCenters= 0;
  for( ; ax<nAxes && !dataFile.eof(); ++ax )
  {
    dataFile>>sizes[ax];
    nCenters+= sizes[ax];
  }
  if( ax<nAxes ) return 0;

  // extract bin centers
  ax= 0;
  double *fBinCenters= new double[nCenters];
  double *binCenters= fBinCenters;
  for( ; ax<nAxes && !dataFile.eof(); ++ax )
  {
    int bin= 0;
    for( ; bin<sizes[ax] && !dataFile.eof(); ++bin, ++binCenters )
      dataFile>>binCenters[0];
    if( bin!=sizes[ax] ) break;
  }
  if( ax==nAxes ) return new PDF_Axes(nAxes, sizes, fBinCenters);
  delete sizes;
  delete fBinCenters;
  return 0;
}

PDF_Axes::~PDF_Axes()
{
    if( m_BinCenters ) delete []m_BinCenters;
    if( m_Sizes )      delete []m_Sizes;
    if( m_Bins )       delete []m_Bins;
}

double PDF_Axes::getBinCenter(int ax, int idAx) const
{
  const double *axData= m_BinCenters;
  for( int i= 0; i<ax;  axData+= m_Sizes[i++] );
  if(abs(idAx)>m_Sizes[ax]) return 0;
  if(idAx<0) return axData[m_Sizes[ax]+idAx];
  return axData[idAx];
}

const double* PDF_Axes::getBinCenters(int ax, bool bin) const
{
  const double *axData= m_BinCenters;
  for( int i= 0; i<ax;  axData+= m_Sizes[i++] );
  if(bin) axData+= m_Bins[ax];
  return axData;
}

bool PDF_Axes::findBin(int axis, double data) const
{
  // finds the bin corresponding to the double  on the axis int
  // returns true if not found: point outside our phase space
  const double *axData= getBinCenters(axis, false);
  int nabove= m_Sizes[axis]+1;
  int nbelow= 0;
  int middle= 0;

  while(nabove-nbelow>1) 
  {
    middle = (nabove+nbelow)/2;
    if (data==axData[middle-1])
    {
      nbelow= middle;
      break;
    }
    if (data<axData[middle-1]) nabove = middle;
    else nbelow = middle;
  }

  m_Bins[axis]= --nbelow;
  // check if data is not inside axis bounds
  if( nbelow<0 || nbelow>=m_Sizes[axis]-1 ) return true;
  return false;
}

bool PDF_Axes::findBin(double* data, int &result) const {
  // finds the pahse space bin corresponding to the parameter vector double*
  // returns true if not found: point outside our phase space
  result= 0;
  int width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax ){
    if( findBin( ax, data[ax] ) ) return true;
    result+= width*m_Bins[ax];
    width*= m_Sizes[ax];
  }
  return false;
}
    

bool PDF_Axes::init(const double *data)
{
  for( int ax= 1; ax<m_Naxes; ++ax )
    if( findBin( ax, data[ax] ) ) return true;
  return false;
}

int PDF_Axes::calculateNeigbouringBin(int n) const
{
  int id= 0, width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax )
  {
    id+= width*(m_Bins[ax]+((n&(1<<ax))>0));
    width*= m_Sizes[ax];
  }
  return id;
}

int PDF_Axes::getNbins(void) const
{
  int width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax )
    width*= m_Sizes[ax];
  return width;
}

PDF_Data *PDF_Data::read(std::ifstream &dataFile,
                         int npar,
                         const PDF_Axes*axes)
{
  if( !axes ) return 0;

  int nBins= axes->getNbins()*npar;
  float *fData= new float[nBins];
  float *data= fData, *end= fData+nBins;
  for( ; data!=end && !dataFile.eof(); ++data ) dataFile>>*data;

  if( data==end ) 
    return new PDF_Data(fData, npar, axes);
  else
  { 
    delete fData; 
    return 0; 
  }
}

PDF_Data::~PDF_Data()
{ 
    m_Axes = 0;
    if (m_Data) delete []m_Data;
}

bool PDF_Data::evaluateParameters(const double *data, double *par ) const 
{
    if( m_Axes->findBin(0, data[0]) ) return true;

    // get weights
    double *weights= new double[m_Axes->getNaxes()];
    for( int ax= m_Axes->getNaxes()-1; ax>=0; --ax )
    {
        const double* binCenter = m_Axes->getBinCenters(ax, true);
        weights[ax] = (data[ax] - binCenter[0]) / (binCenter[1] - binCenter[0]);
    }
  
    // linear interpolation
    float* values = new float[m_Npar];
    for( int ax= 0; ax < m_Npar; ++ax ) par[ax]= 0.;
    for( int n = (1 << m_Axes->getNaxes()) - 1; n >= 0; --n )
    {
        // find position for neighbouring bin
        // and multiply values by the corresponding weight 
        // code is: n&(1<<ax) means neighbour at +1 for axis number ax
        memcpy(values, m_Data + m_Axes->calculateNeigbouringBin(n) * m_Npar, 4 * m_Npar);
        for( int id = m_Axes->getNaxes()-1; id >= 0; --id ) 
        {
            double w = (1 << id) & n ? weights[id] : 1 - weights[id];
            for( int iPar = 0; iPar < m_Npar; ++iPar ) values[iPar] *= w;
        }
        for( int iPar= 0; iPar < m_Npar; ++iPar ) par[iPar] += values[iPar];
    }

    delete []weights;
    delete []values;
    return false;
}

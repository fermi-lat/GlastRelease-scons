#include "LikelihoodTool.h"
#include "facilities/Util.h"
#include <fstream>

StatusCode LikelihoodTool::readPDFparameters( MsgStream &log,
                                       std::string &dataFileName ){
    // This function does following initialization actions:
    //    - extracts PDF parameters from txt file
    StatusCode sc = StatusCode::SUCCESS;
	  log << MSG::DEBUG <<"Initializing LikelihoodTool with file "<<dataFileName<<endreq;

    // Read in the parameters from the data file
    facilities::Util::expandEnvVar(&dataFileName);
    std::ifstream dataFile(dataFileName.data());
    if( !dataFile ){
      log<<MSG::ERROR<<" Unable to open data file: "<<dataFileName<<endreq;
      return StatusCode::FAILURE;
    }
    
    int nPar;
    std::string line;
    getline(dataFile, line);

    m_PDFAxes= 0;
    m_PDFCol= 0;
    int pdf= 0;
    int errflag;
    while( !dataFile.eof() ){
      getline(dataFile, line);
      if( (line=="") || (line=="\n") ) continue;
      else if( line.size()>6 && (line.substr(0, 7)=="COMMENT") ) continue; 
      else if( line.size()>6 && (line.substr(0, 6)=="#PDFs:") ){
        m_PDFCol= new PDF_Data*[m_Npdf];
        for( int pdfs= 0; pdfs<m_Npdf; ++pdfs ) m_PDFCol[pdfs]= 0;
        sscanf(line.data(), "#PDFs: %d #Parameters: %d", &m_Npdf, &nPar);
        log<<MSG::DEBUG<<"File must contain "<<m_Npdf<<", for "<<nPar
                       <<" parameters"<<endreq;
        if( nPar>0 && m_Npdf>0 ) continue;
      } else if( line=="AXES:" ){
        if( !m_PDFAxes ){
          log<<MSG::DEBUG<<"Extracting Axes"<<endreq;
          m_PDFAxes= new PDF_Axes(log, dataFile, errflag);
          if( !errflag ) continue;
        }
      } else if( line.size()>4 && line.substr(0,4)=="PDF:" ){
        if( !m_PDFCol ) {
          log<<MSG::ERROR<<"Line \"#PDFs: <number of PDFs> "
                         <<"Parameters: <number of parameters for PDF function>"
                         <<"\" must come  before 1st PDF";
        } else if( pdf<m_Npdf ) {
          log<<MSG::DEBUG<<"Extracting PDF Parameters for PDF #"
                         <<pdf<<endreq;
          m_PDFCol[pdf++]= 
            new PDF_Data(log, dataFile, m_PDFAxes, nPar, errflag);
          if( !errflag ) continue;
        }
      }
      log<<MSG::ERROR<<"Incorrect data in file: "<<dataFileName<<endreq
                     <<"Stopping on data line: \""<<line<<"\""<<endreq;
      dataFile.close();
      return StatusCode::FAILURE;
    }
    dataFile.close();
    if( pdf<m_Npdf || !m_PDFAxes ){
      log<<MSG::ERROR<<"Missing data in file: "<<dataFileName<<endreq;
      return StatusCode::FAILURE;
    }
    return sc;
}


StatusCode LikelihoodTool::calculateEvent( int iEvtPDF,
                              double methodMinEnergy, double methodMaxEnergy,
                              double *pdfVariables, double *pdfDataPoint,
                              double &recEnergy, double &recEnergyWidth )
//Purpose and method:
//   This function performs:
//        Evaluation of the energy for a set of parameters (pdfDataPoint),
//      specifying the form of a function (pdfFCN) which is evaluated at a
//      given point (pdfVariables, reconstruction energy). 
//        The function parameters are evaluated using linear in interpolation
//      of data points using eventPDF. 
//        The energy is estimated as the PDF's MPV for that given set of data
//      and variables. The quality of the reconstruction is estimated as the
//      FWHM of the PDF.
//        Both values are estimated recusively, calculating 10 points and
//      choosing a range around that turn's best estimate.
//      
//   The main actions are:
//      - estimate particle energy by TKR correction method (CALCULATE)
//      - estimate FWHM of the PDF for set a of parameters: 
//          pdfDataPoint && reconstructerd energy            (QUALITY)
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
  StatusCode sc= StatusCode::SUCCESS;
  const PDF_Data *eventPDF= m_PDFCol[iEvtPDF];
  double minE= pdfVariables[0]*.5;
  double maxE= pdfVariables[0]*5.;
  if( minE<methodMinEnergy ) minE= methodMinEnergy;
  if( maxE>methodMaxEnergy ) maxE= methodMaxEnergy;

  // next variables are for the estimation of the quality of the 
  // reconstruction: m_widthEnergyCorr
  double fwhmLimits[2]= { minE, maxE };
  double fwhmBinWidth[2];
  fwhmBinWidth[0]= fwhmBinWidth[1]= maxE-minE;

  recEnergy= minE;
  // stores reconstructed energy pdf value
  double maxPdf= -1.;

  double *pdfParameters= new double[eventPDF->getNpar()];
  pdfDataPoint[0]= minE;
  m_PDFAxes->init(pdfDataPoint);
  for( double bW= (maxE-minE)/10.; bW>.1; bW= (maxE-minE)/10. ){
    int errCalls= 0;
    for( pdfDataPoint[0]= minE; pdfDataPoint[0]<maxE; pdfDataPoint[0]+= bW ){
      // get log normal parameters
      if( eventPDF->evaluateParameters(pdfDataPoint, pdfParameters) ) 
      { ++errCalls; continue; }

      // calculate pdf value with these parameters
      double pdfVal= pdfFCN(pdfVariables, pdfParameters);
      if( pdfVal>maxPdf ) { maxPdf= pdfVal; recEnergy= pdfDataPoint[0]; } 
      else if( pdfVal<maxPdf*.5 ){
        int iFWHM= recEnergy<pdfDataPoint[0];
        if(   (iFWHM ^ (pdfDataPoint[0]>fwhmLimits[iFWHM]))
           || (iFWHM ^ (recEnergy<fwhmLimits[iFWHM])) ) {
          fwhmLimits[iFWHM]= pdfDataPoint[0];
          fwhmBinWidth[iFWHM]= bW;
        }
      }
    }
    if( errCalls==10 ) {
      recEnergy= (double) cutNOPARAMETERS;
      recEnergyWidth= 1.;
      return sc;
    }

    minE= recEnergy-bW;
    maxE= recEnergy+bW;
    if( minE<methodMinEnergy ) minE= methodMinEnergy;
    if( maxE>methodMaxEnergy ) maxE= methodMinEnergy;
  }

  // QUALITY
  // now looking for FWHM
  // fwhmLimits is an outer limit for an energy such that pdfFCN(x)<.5*maxPdf
  // fwhmBinWidth gives us the inner limit


  recEnergyWidth= 0.;
  for( int iFWHM= 0; iFWHM<2; ++iFWHM ){
    // reached limits of reconstruction range
    if( fwhmLimits[iFWHM]==0. ) {
      recEnergyWidth+= 1.;
      continue;
    }

    for( double bW= fwhmBinWidth[iFWHM]/10.; bW>.01; bW= (maxE-minE)/10. ){
      if( iFWHM ) { maxE= fwhmLimits[1]; minE= maxE-bW; }
      else { minE= fwhmLimits[0]; maxE= minE+bW; }

      int errCalls= 0;
      for( pdfDataPoint[0]= minE; pdfDataPoint[0]<maxE; pdfDataPoint[0]+= bW ){
        if( eventPDF->evaluateParameters(pdfDataPoint, pdfParameters) ) { 
          ++errCalls; 
          continue; 
        } 
        if(    (pdfFCN(pdfVariables, pdfParameters)<maxPdf*.5)
            && (iFWHM ^ (pdfDataPoint[0]>fwhmLimits[iFWHM])) )
          fwhmLimits[iFWHM]= pdfDataPoint[0];
      }
      if( errCalls==10 ){ 
        fwhmLimits[iFWHM]= 0.;
        break;
      }
    }
    recEnergyWidth+= .5*fabs(recEnergy-fwhmLimits[iFWHM])/recEnergy;
  }

  delete []pdfParameters;
  return sc;
}

double LikelihoodTool::pdfFCN( double x[2], double parameters[5] ) const{
  // calculate LogNormal(x[1]*parameters[0]+x[0]):
  // \f[ LogNormal(x, parameters(recEnergy)) \f]
  // if( par[3]==0. ) return 0.; that shouldn't happen
  // \f[ LogNormal(x)= N\,\exp(-\frac{1}{2}(\frac{\ln(1+\tau (x-\mu) \frac{\sinh({\tau}\sqrt{\ln\,4})}{ 2.36\beta{\tau}\,\sqrt{\ln\,4}})}\tau)^2+\tau^2) \f]
  double result= (x[0]+x[1]*parameters[0]-parameters[2])/parameters[3];
  if( parameters[4]!=0. ){
    double val= 1.17741002251547466*parameters[4]; //sqrt(log(4.))*\tau
    result*= parameters[4]*(sinh(val)/val);
    result= log(1+result)/parameters[4]; 
  }
  result*= -result*.5;
  result= parameters[1]*exp(result+parameters[4]*parameters[4]);
  return isnan(result)?0.:result;
}

StatusCode LikelihoodTool::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  if( m_PDFCol ){
    for( int plane= 0; plane<m_Npdf; ++plane )
      if( m_PDFCol[plane] ){ delete m_PDFCol[plane]; m_PDFCol[plane]= 0; }
    delete []m_PDFCol; 
  }
  if( m_PDFAxes ) delete m_PDFAxes;
  return sc;
}

PDF_Axes::PDF_Axes(MsgStream &log, std::ifstream &dataFile, int &ierr){
  m_Sizes= 0;
  m_BinCenters= 0;
  ierr= true;

  dataFile>>m_Naxes;
  m_Sizes= new int[m_Naxes];
  m_Bins= new int[m_Naxes];

  // next are axis sizes
  int ax= 0;
  int nCenters= 0;
  log<<MSG::DEBUG<<"Reading Axes ";
  for( ; ax<m_Naxes && !dataFile.eof(); ++ax ){
    dataFile>>m_Sizes[ax];
    nCenters+= m_Sizes[ax];
    log<<MSG::DEBUG<<"\t"<<m_Sizes[ax];
  }
  log<<MSG::DEBUG<<endreq;
  log<<MSG::DEBUG<<m_Naxes<<" Axes for a total of "<<getNbins()<<" bins"
                 <<endreq;
  if( ax<m_Naxes ) return;
  m_BinCenters= new double[nCenters];
  

  // next are bin centers
  double *binCenters= m_BinCenters;
  ax= 0;
  for( ; ax<m_Naxes && !dataFile.eof(); ++ax ){
    int bin= 0;
    log<<MSG::DEBUG<<"Reading Centers for Axis "<<ax<<":"<<endreq;
    for( ; bin<m_Sizes[ax] && !dataFile.eof(); ++bin, ++binCenters ){
      dataFile>>binCenters[0];
      log<<MSG::DEBUG<<"("<<bin<<":"<<binCenters[0]<<")";
    }
    log<<MSG::DEBUG<<endreq;
    if( bin!=m_Sizes[ax] ){
      log<<MSG::ERROR<<"Wrong Number of values for axis"<<ax<<endreq;
      return;
    }
  }
  log<<MSG::DEBUG<<"Read "<<(binCenters-m_BinCenters)<<" values"<<endreq;

  ierr= ax<m_Naxes;
}

PDF_Axes::~PDF_Axes(){
  if( m_BinCenters ) delete []m_BinCenters;
  if( m_Sizes ) delete []m_Sizes;
  if( m_Bins ) delete []m_Bins;
  m_Sizes= 0;
  m_BinCenters= 0;
}

const double* PDF_Axes::getBinCenters(int ax, bool bin) const {
  const double *axData= m_BinCenters;
  for( int i= 0; i<ax;  axData+= m_Sizes[i++] );
  if(bin) axData+= m_Bins[ax];
  return axData;
}

bool PDF_Axes::findBin(int axis, double data) const {
  const double *axData= getBinCenters(axis, false);
  int nabove= m_Sizes[axis]+1;
  int nbelow= 0, middle= 0;
  while(nabove-nbelow>1) {
    middle = (nabove+nbelow)/2;
    if (data==axData[middle-1]) { nbelow= middle; break; }
    if (data<axData[middle-1]) nabove = middle;
    else nbelow = middle;
  }
  m_Bins[axis]= --nbelow;
  // check if data is not inside axis bounds
  if( nbelow<0 || nbelow>=m_Sizes[axis]-1 ) return true;
  return false;
}

bool PDF_Axes::init(const double *data) {
  for( int ax= 1; ax<m_Naxes; ++ax )
    if( findBin( ax, data[ax] ) ) return true;
  return false;
}

int PDF_Axes::calculateNeigbouringBin(int n) const {
  int id= 0, width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax ){
    id+= width*(m_Bins[ax]+(n&(1<<ax)));
    width*= m_Sizes[ax];
  }
  return id;
}

int PDF_Axes::getNbins(void) const {
  int width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax ) width*= m_Sizes[ax];
  return width;
}

PDF_Data::PDF_Data(MsgStream &log, std::ifstream &dataFile, PDF_Axes*axes,
                   int npar, int &ierr):
    m_Npar(npar), m_Axes(axes){
  ierr= true;
  if( !m_Axes ) return;
  int nBins= m_Axes->getNbins()*m_Npar;
  log<<MSG::DEBUG<<"Reading "<<m_Npar<<" parameters per bin for "
                 <<nBins/m_Npar<<" bins"<<endreq;
  m_Data= new float[nBins];
  float *data= m_Data, *end= m_Data+nBins;
  for( ; data!=end && !dataFile.eof(); ++data ) dataFile>>*data;
  log<<MSG::DEBUG<<"Read "<<int(data-m_Data)/m_Npar<<" bins"<<endreq;
  ierr= data!=end;
}

PDF_Data::~PDF_Data(){ 
  m_Axes=0;
  if( m_Data ) delete []m_Data;
}

bool PDF_Data::evaluateParameters(const double *data, double *par ) const {
  if( m_Axes->findBin(0, data[0]) ) return true;

  // get weights
  double *weights= new double[m_Axes->getNaxes()];
  for( int ax= m_Axes->getNaxes()-1; ax>=0; --ax ){
    const double *binCenter= m_Axes->getBinCenters(ax, true);
    weights[ax]= (data[ax]-binCenter[0])/(binCenter[1]-binCenter[0]);
  }
  
  // linear interpolation
  float *values= new float[m_Npar];
  for( int ax= 0; ax<m_Npar; ++ax ) par[ax]= 0.;
  for( int n= (1<<m_Axes->getNaxes())-1; n>=0; --n ){
    // find position for neighbouring bin
    // and multiply values by the corresponding weight 
    // code is: n&(1<<ax) means neighbour at +1 for axis number ax
    memcpy(values, m_Data+m_Axes->calculateNeigbouringBin(n)*m_Npar, 4*m_Npar);
    for( int id= m_Axes->getNaxes()-1; id>=0; --id ) {
      double w= (1<<id)&n? weights[id]:1-weights[id];
      for( int iPar= 0; iPar<m_Npar; ++iPar ) values[iPar]*= w;
    }
    for( int iPar= 0; iPar<m_Npar; ++iPar ) par[iPar]+= values[iPar];
  }

  delete []weights;
  delete []values;
  return false;
}

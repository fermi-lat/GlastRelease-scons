#ifndef __CalLikelihoodPDF__CXX
#include "CalLikelihoodPDFParameters.h"
#include "GaudiKernel/MsgStream.h"
#include <fstream>

/******************************************************************************/
/***********************  PDFParameters ********************************************/
/******************************************************************************/
PDFParameters::PDFParameters(std::ifstream &dataFile, PDFGrid*axes, int nPar, 
                   MsgStream &log)
  :m_Npar(nPar), m_Data(0), m_Grid(axes)
{
  if( !axes  ){
    log<<MSG::ERROR<<"ERROR Creating PDFParameters: Axes is NULL"<<std::endl;
    return;
  }
  int nBins= axes->getNbins()*getNpar();
  double *fData= new double[nBins];
  double *data= fData, *end= fData+nBins;
  for( ; data!=end && !dataFile.eof(); ++data ) dataFile>>*data;

  if( data!=end ){
    log<<MSG::ERROR<<"ERROR Creating PDFParameters: Missing Data"<<std::endl;
    delete []fData;
  } else m_Data= fData;
}

PDFParameters::~PDFParameters()
{ 
    m_Grid= 0;
    if (m_Data) delete []m_Data;
}

bool PDFParameters::interpolation(const double data[], double par[], bool ini) {
    if( ini ) m_Grid->initialise(data);
    if( m_Grid->findBin(0, data[0]) ) return true;
    // get weights
    double *weights= new double[m_Grid->getNaxes()];
    for( int ax= m_Grid->getNaxes()-1; ax>=0; --ax )
    {
        const double* binCenter = m_Grid->getBinCenters(ax, true);
        weights[ax] = (data[ax]-binCenter[0]) / (binCenter[1]-binCenter[0]);
    }
  
    // linear interpolation
    double* values = new double[getNpar()];
    for( int ax= 0; ax < getNpar(); ++ax ) par[ax]= 0.;
    for( int n = (1 << m_Grid->getNaxes()) - 1; n >= 0; --n )
    {
        // find position for neighbouring bin
        // and multiply values by the corresponding weight 
        // code is: n&(1<<ax) means neighbour at +1 for axis number ax
        memcpy(values, m_Data + m_Grid->calculateNeigbouringBin(n)*getNpar(),
               sizeof(double)*getNpar());
        for( int id = m_Grid->getNaxes()-1; id >= 0; --id ) 
        {
            double w = (1 << id) & n ? weights[id] : 1 - weights[id];
            for( int iPar = 0; iPar < getNpar(); ++iPar ) values[iPar] *= w;
        }
        for( int iPar= 0; iPar < getNpar(); ++iPar ) par[iPar] += values[iPar];
    }

    delete []weights;
    delete []values;
    return false;
}

/******************************************************************************/
/***********************  PDFGrid *********************************************/
/******************************************************************************/
PDFGrid::PDFGrid(std::ifstream &dataFile, MsgStream &log)
  :m_BinCenters(0), m_Sizes(0), m_Bins(0), m_Naxes(0)
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
  if( ax<nAxes ) {
    log<<MSG::ERROR<<"ERROR Creating PDFGrid: Missing Data"<<std::endl;
    return;
  }

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
  if( ax==nAxes ){
    m_Naxes= nAxes;
    m_BinCenters= fBinCenters;
    m_Sizes= sizes;
    m_Bins= new int[nAxes];
    return;
  } else {
    log<<MSG::ERROR<<"ERROR Creating PDFGrid: Missing Data"<<std::endl;
    return;
  }
  delete []sizes;
  delete []fBinCenters;
  return;
}

PDFGrid::~PDFGrid()
{
    if( m_BinCenters ) delete []m_BinCenters;
    if( m_Sizes )      delete []m_Sizes;
    if( m_Bins )       delete []m_Bins;
}

const double* PDFGrid::getBinCenters(int ax, bool bin) const
{
  const double *axData= m_BinCenters;
  for( int i= 0; i<ax;  axData+= m_Sizes[i++] );
  if(bin) axData+= m_Bins[ax];
  return axData;
}

double PDFGrid::getBinCenter(int ax, int idAx) const
{
  const double *axData= m_BinCenters;
  for( int i= 0; i<ax;  axData+= m_Sizes[i++] );
  if(abs(idAx)>m_Sizes[ax]) return 0;
  if(idAx<0) return axData[m_Sizes[ax]+idAx];
  return axData[idAx];
}

bool PDFGrid::findBin(int axis, double data) const
{
  const double *axData= getBinCenters(axis, false);
  int nabove= m_Sizes[axis]+1;
  int nbelow= 0;

  while(nabove-nbelow>1) 
  {
    int middle = (nabove+nbelow)/2;
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

bool PDFGrid::findBin(double* data, int &result) const {
  result= 0;
  int width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax ){
    if( findBin( ax, data[ax] ) ) return true;
    result+= width*m_Bins[ax];
    width*= m_Sizes[ax];
  }
  return false;
}

bool PDFGrid::initialise(const double data[])
{
  for( int ax= 1; ax<m_Naxes; ++ax )
    if( findBin( ax, data[ax] ) ) return true;
  return false;
}

int PDFGrid::calculateNeigbouringBin(int n) const
{
  int id= 0, width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax )
  {
    id+= width*(m_Bins[ax]+((n&(1<<ax))>0));
    width*= m_Sizes[ax];
  }
  return id;
}

int PDFGrid::getNbins(void) const
{
  int width= 1;
  for( int ax= m_Naxes-1; ax>=0; --ax )
    width*= m_Sizes[ax];
  return width;
}
#endif

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Digi/TkrDigi.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "facilities/Util.h"

#include <CalRecon/ICalEnergyCorr.h>

#include <fstream>

class PDF_Axes 
{
public:
    PDF_Axes(MsgStream&, std::ifstream&, int&);
   ~PDF_Axes();

    bool init( const double* );
    int getNbins(void) const;
    int getNaxes(void) const { return m_Naxes; }
    int getSize(int ax) const { return m_Sizes[ax]; }
    const double* getBinCenters(int, bool) const;

    bool findBin(int, double) const;
    int calculateNeigbouringBin(int) const;
private:
    double *m_BinCenters;
    int *m_Sizes;
    int *m_Bins;
    
    int m_Naxes;
};

class PDF_Data 
{
public:
    PDF_Data(MsgStream&, std::ifstream&, PDF_Axes*, int, int&);
   ~PDF_Data();

    int getNpar( void ) const { return m_Npar; }
    bool evaluateParameters(const double*, double*) const;
private:
    float *m_Data;
    int m_Npar;
    const PDF_Axes *m_Axes;
};

/**   
* @class CalTkrLikelihoodTool
*
* Algorithm for correction of energy degradation in the tracker by correlating
* the energy in the CAL with the number of hit TKR strips.  
*
*
*/

class CalTkrLikelihoodTool : public AlgTool, virtual public ICalEnergyCorr
{
public:
    typedef enum { cutNOTKRREC     = -1, 
                   cutMINCALENERGY = -2, 
                   cutMAXCALENERGY = -3,
                   cutSLOPE        = -4, 
                   cutVERTEX       = -5, 
                   cutPOSITION     = -6,
                   cutNOPARAMETERS = -7 }
                 LikelihoodTool_Cuts;

    //! constructor
    CalTkrLikelihoodTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent);
    virtual ~CalTkrLikelihoodTool(){}
    
    StatusCode initialize();

    StatusCode finalize();

    //! Tracker energy degradation correction using number of TKR hit strips
    /*! This method uses the correlation between the energy \"lost\" in the 
    * tracker and the energy deposited in the calorimeter.
    * We used the Monte Carlo simulation of the LAT to determine this correlation
    * at several energies, from 50 MeV up to 1 GeV, and angles from 0 to 32\deg. 
    * For one particular incident energy and angle, the bidimensionnal
    * distribution of  the  number of hit strips and the energy deposited in the 
    * CAL can be characterised by the 1D distribution:
    * \f[ E_{CAL} + \alpha TkrHits \f]
    * where  \f$\alpha\f$ is been optimised so as to obtain the narrowest such
    * distribution, normalised to a probability and with its MPV at the incident
    * energy.
    * These distributions can be used to defined a probability density function.
    * The reconstructed energy for a given event then becomes the one maximising
    * the probability, for a reconstruced direction, CAL raw energy,...
    *
    * \par The method takes 4 arguments:
    * \param eTotal Total energy measured in the calorimeter in MeV
    * \param nHits  Total number of hit strips in the CAL
    * \param vertex[3] reconstructed vertex position  
    * \param dir[3] reconstructed direction
    *
    *\return Corrected energy in MeV
    *
    *\warning needs TKR reconstruction
    *
    *\author
    */
           
    Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrVertex* );
    
private:
    // basic building brick for any PDF
    StatusCode calculateEvent( int, double, double, double*,
                               double*, double&, double& );
    StatusCode readPDFparameters( MsgStream &lm, std::string& );
    double integratedDistance2TowerSide(double, double) const;
    double pdfFCN(double[2], double[5]) const;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;

    /// Detector Service
    IGlastDetSvc*     m_detSvc; 

    // data tables
    PDF_Axes *m_PDFAxes;
    PDF_Data **m_PDFCol;
    int m_Npdf;

    // CAL origin, tower pitch, and gap between a CAL and its tower,
    // needed for the cuts
    double m_zOriginCAL;
    double m_pitchTOWER;
    double m_halfwidthCAL;
    double m_heightCAL;
    
    /// input data file containing log normal parameters
    std::string m_dataFile;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalTkrLikelihoodTool) ;

CalTkrLikelihoodTool::CalTkrLikelihoodTool( const std::string& type, 
                                            const std::string& name, 
                                            const IInterface* parent)
                                          : AlgTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ICalEnergyCorr>(this);
    declareProperty("dataFile",
                    m_dataFile="$(CALRECONROOT)/xml/CalTkrLikelihood.data");
};

StatusCode CalTkrLikelihoodTool::initialize()
{
    // This function does following initialization actions:
    //    - extracts geometry constants from xml file using GlastDetSvc
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
      log << MSG::INFO << "Initializing CalTkrLikelihoodTool" <<endreq;

    typedef std::map<double*,std::string> PARAMAP;
    PARAMAP param;
    param[&m_pitchTOWER]  = std::string("towerPitch");
    param[&m_halfwidthCAL]= std::string("CALModuleWidth");
    param[&m_heightCAL]   = std::string("CALModuleHeight");

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

    m_zOriginCAL    = -47.395;
    m_halfwidthCAL *= .5;
    
    // Read in the parameters from the data file
    readPDFparameters( log, m_dataFile );

    return sc;
}

StatusCode CalTkrLikelihoodTool::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    if( m_PDFCol )
    {
        for( int plane= 0; plane<m_Npdf; ++plane )
        {
            if( m_PDFCol[plane] ){ delete m_PDFCol[plane]; m_PDFCol[plane]= 0; }
        }
        delete []m_PDFCol; 
    }
  
    if( m_PDFAxes ) delete m_PDFAxes;
  
    return sc;
}

Event::CalCorToolResult* CalTkrLikelihoodTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrVertex* vertex)
//Purpose and method:
//
//   This function performs:
//   The main actions are:
//      - check wheter the event meets basic requirements (CUTS)
//      - calculate energy by TKR correction method using LikelihoodTool 
// 
// TDS input: CalCluster (and TkrDigi?)
// TDS output: CalClusters
{
    Event::CalCorToolResult* corResult = 0;

    // if reconstructed tracker data doesn't exist or number of tracks is 0:
    if (vertex == 0)
    {
        //cluster->setEnergyCalTkrLikelihood((double) LikelihoodTool::cutNOTKRREC, 1.);
        return corResult;
    }

    const Vector& trackDirection = vertex->getDirection();
    const Point&  trackPosition  = vertex->getPosition();
      
    // CUTS
    // this checks whether a set of PDF parameters exist for this event's
    // energy, direction and point of impact.
    // if not, the setEnergyCorr value corresponds to the
    // defined CalTkrLikelihoodToolFlags values
  
    // Energy:
    double pdfVariables[2] = {cluster->getCalParams().getEnergy(), 0.};
    double pdfDataPoint[2] = {0., fabs(trackDirection.z())};

    if( pdfVariables[0] < 20. ) 
    {
        //cluster->setEnergyCalTkrLikelihood((double) LikelihoodTool::cutMAXCALENERGY, 1.);
        return corResult;
    }

    if( pdfVariables[0] > 1900. ) 
    {
        //cluster->setEnergyCalTkrLikelihood((double) LikelihoodTool::cutMAXCALENERGY, 1.);
        return corResult;
    }
  
    // direction: slope must be above \f$cos(32\circ)$\f
    if( fabs(trackDirection.z()) < .8480481 )
    { 
        //cluster->setEnergyCalTkrLikelihood((double) LikelihoodTool::cutSLOPE, 1.);
        return corResult; 
    }

    // z-position:
    int vertexPos =  int((trackPosition[2]-108.)*3.2e-2);
    if( vertexPos < 0 || vertexPos > 15 ) 
    { 
        //cluster->setEnergyCalTkrLikelihood((double) LikelihoodTool::cutVERTEX, 1.);
        return corResult; 
    }

    // x,y-position: 1. find the track position on top of the CAL
    //               2. calculate its distance to parallel sides of the CAL,
    //    integrated along the track, normalized to 1.
    //               3. use the geometric mean of values on X and Y axis as a
    //    basis for geometric cut.
    double geometricCut = 1.;
    for( int ax= 0; ax<2; ++ax )
    {
        double slope = -trackDirection[ax] / trackDirection[2];
        double posCAL = trackPosition[ax];

        posCAL       += (trackPosition[2] - m_zOriginCAL) * slope;
        geometricCut *= integratedDistance2TowerSide(posCAL, slope);
    }
    geometricCut= sqrt(geometricCut);

    // correcting for event slope
    double thetaNorm = .707106781186547462 * sqrt(1. / (trackDirection[2] * trackDirection[2]) - 1.);

    thetaNorm    *= m_heightCAL / (2 * m_halfwidthCAL); 
    geometricCut /= (1 - thetaNorm);
  
    if( geometricCut<.01 )
    {
        //cluster->setEnergyCalTkrLikelihood((double) cutPOSITION, 1.);
        return corResult;
    }
    int geometricCutindex= (geometricCut > .35) + (geometricCut > .75);

    // CALCULATE
    //    - get number of hits in TKR
    Event::TkrDigiCol* tkrDigiData = 
                            SmartDataPtr<Event::TkrDigiCol>(m_dataSvc, EventModel::Digi::TkrDigiCol); 
    for (Event::TkrDigiCol::iterator it = tkrDigiData->begin(); it != tkrDigiData->end(); ++it)
         pdfVariables[0] += (double)(*it)->getNumHits();

    double recEnergy      = 0.;
    double recEnergyWidth = 1.;

    StatusCode sc = calculateEvent( vertexPos + geometricCutindex * 15,
                                   50., 2000., pdfVariables, pdfDataPoint,
                                   recEnergy, recEnergyWidth );

    if (sc.isSuccess())
    {
        // Create a CalCorToolResult object to hold the information
        corResult = new Event::CalCorToolResult();

        corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
        corResult->setCorrectionName(type());
        corResult->setParams(cluster->getCalParams());
        corResult->setChiSquare(1.);
        corResult->insert(Event::CalCorEneValuePair("RecEnergy",recEnergy));
        corResult->insert(Event::CalCorEneValuePair("RecEnergyWidth",recEnergyWidth));
    }

    return corResult;
}

StatusCode CalTkrLikelihoodTool::calculateEvent(int    iEvtPDF,
                                                double methodMinEnergy, 
                                                double methodMaxEnergy,
                                                double *pdfVariables, 
                                                double *pdfDataPoint,
                                                double &recEnergy, 
                                                double &recEnergyWidth )
{
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

    const PDF_Data* eventPDF = m_PDFCol[iEvtPDF];
    double          minE     = pdfVariables[0] * .5;
    double          maxE     = pdfVariables[0] * 5.;

    if (minE < methodMinEnergy ) minE = methodMinEnergy;
    if (maxE > methodMaxEnergy ) maxE = methodMaxEnergy;

    // next variables are for the estimation of the quality of the 
    // reconstruction: m_widthEnergyCorr
    double fwhmLimits[2] = { minE, maxE };
    double fwhmBinWidth[2];

    fwhmBinWidth[0] = fwhmBinWidth[1] = maxE - minE;

    recEnergy = minE;
    // stores reconstructed energy pdf value
    double maxPdf = -1.;

    double* pdfParameters = new double[eventPDF->getNpar()];
    pdfDataPoint[0] = minE;
    m_PDFAxes->init(pdfDataPoint);

    for(double bW = (maxE - minE) / 10.; bW > .1; bW = (maxE - minE) / 10. )
    {
        int errCalls = 0;
        for(pdfDataPoint[0] = minE; pdfDataPoint[0] < maxE; pdfDataPoint[0] += bW)
        {
            // get log normal parameters
            if( eventPDF->evaluateParameters(pdfDataPoint, pdfParameters) ) 
            { 
                ++errCalls; 
                continue;
            }

            // calculate pdf value with these parameters
            double pdfVal= pdfFCN(pdfVariables, pdfParameters);
            if( pdfVal>maxPdf ) 
            { 
                maxPdf    = pdfVal; 
                recEnergy = pdfDataPoint[0]; 
            } 
            else if (pdfVal < maxPdf * .5 )
            {
                int iFWHM= recEnergy<pdfDataPoint[0];
                if(   (iFWHM ^ (pdfDataPoint[0]>fwhmLimits[iFWHM]))
                   || (iFWHM ^ (recEnergy<fwhmLimits[iFWHM])) ) 
                {
                    fwhmLimits[iFWHM]   = pdfDataPoint[0];
                    fwhmBinWidth[iFWHM] = bW;
                }
            }
        }
        if( errCalls==10 ) 
        {
            recEnergy      = (double) cutNOPARAMETERS;
            recEnergyWidth = 1.;
            return sc;
        }

        minE = recEnergy - bW;
        maxE = recEnergy + bW;
        if( minE<methodMinEnergy ) minE = methodMinEnergy;
        if( maxE>methodMaxEnergy ) maxE = methodMinEnergy;
    }

    // QUALITY
    // now looking for FWHM
    // fwhmLimits is an outer limit for an energy such that pdfFCN(x)<.5*maxPdf
    // fwhmBinWidth gives us the inner limit


    recEnergyWidth = 0.;
    for( int iFWHM= 0; iFWHM<2; ++iFWHM )
    {
        // reached limits of reconstruction range
        if( fwhmLimits[iFWHM]==0. ) 
        {
            recEnergyWidth+= 1.;
            continue;
        }

        for( double bW= fwhmBinWidth[iFWHM]/10.; bW>.01; bW= (maxE-minE)/10. )
        {
            if( iFWHM ) 
            { 
                maxE = fwhmLimits[1]; 
                minE = maxE-bW; 
            }
            else 
            { 
                minE = fwhmLimits[0]; 
                maxE = minE + bW; 
            }

            int errCalls= 0;
            for( pdfDataPoint[0]= minE; pdfDataPoint[0]<maxE; pdfDataPoint[0]+= bW )
            {
                if( eventPDF->evaluateParameters(pdfDataPoint, pdfParameters) ) 
                {
                    ++errCalls; 
                    continue; 
                } 
                if(    (pdfFCN(pdfVariables, pdfParameters)<maxPdf*.5)
                    && (iFWHM ^ (pdfDataPoint[0]>fwhmLimits[iFWHM])) )
                        fwhmLimits[iFWHM]= pdfDataPoint[0];
            }

            if( errCalls==10 )
            { 
                fwhmLimits[iFWHM]= 0.;
                break;
            }
        }

        recEnergyWidth += .5 * fabs(recEnergy - fwhmLimits[iFWHM]) / recEnergy;
    }

    delete []pdfParameters;
    return sc;
}

double CalTkrLikelihoodTool::integratedDistance2TowerSide( double x, double slope ) const 
{
    if( slope == 0.) return 1.;

    if( slope >  0. )
    { 
        x     *= -1.; 
        slope *= -1; 
    }

    bool leaving = false;

    if (x < 0.) x += m_pitchTOWER;
    if (x < 0.) x += m_pitchTOWER;    // TU 5/20/05: Why twice?

    if (x > m_pitchTOWER)
    { 
        x -= m_pitchTOWER; 
        leaving= true; 
    }

    x -= 0.5 * m_pitchTOWER;

    if (fabs(x) > m_pitchTOWER*.5 ) return 0.;
    if (x > m_halfwidthCAL )
    {
        if (leaving) return 0.;
        x -= m_pitchTOWER;
    }

    double height = m_heightCAL;

    if (x < -m_halfwidthCAL)
    {
        if ((height -= (m_halfwidthCAL + x) / slope) <= 0.) return 0.;
        x= -m_halfwidthCAL;
    } 

    double result = 0.;
    double norm   = height * m_halfwidthCAL;

    if (x < 0)
    {
        double z = -x / slope;
        if (-z >= height) 
        {
            result = (m_halfwidthCAL + x - .5 * height * slope) * height / norm;
            return result;
        }
        result  = -(2. * m_halfwidthCAL + x) * z * .5;
        height += z;
        x       = 0.;
    } 
    double z = (m_halfwidthCAL - x) / slope;
    if (-z >= height ) result += (m_halfwidthCAL - x + height * slope * .5) * height;
    else  result -= (m_halfwidthCAL - x) * z * .5;
    result /= norm;
    return result;
}

StatusCode CalTkrLikelihoodTool::readPDFparameters( MsgStream &log, std::string &dataFileName )
{
    // This function does following initialization actions:
    //    - extracts PDF parameters from txt file
    StatusCode sc = StatusCode::SUCCESS;
    
    log << MSG::DEBUG <<"Initializing LikelihoodTool with file "<<dataFileName<<endreq;

    // Read in the parameters from the data file
    facilities::Util::expandEnvVar(&dataFileName);
    std::ifstream dataFile(dataFileName.data());
    if( !dataFile )
    {
      log<<MSG::ERROR<<" Unable to open data file: "<<dataFileName<<endreq;
      return StatusCode::FAILURE;
    }
    
    int nPar;
    std::string line;
    getline(dataFile, line);

    m_PDFAxes = 0;
    m_PDFCol  = 0;
    int pdf = 0;
    int errflag;
    while(!dataFile.eof() )
    {
        getline(dataFile, line);
        if      ((line == "") || (line == "\n"))                        continue;
        else if ( line.size() > 6 && (line.substr(0, 7) == "COMMENT") ) continue; 
        else if ( line.size() > 6 && (line.substr(0, 6) == "#PDFs:") )
        {
            sscanf(line.data(), "#PDFs: %d #Parameters: %d", &m_Npdf, &nPar);
            log << MSG::DEBUG << "File must contain " << m_Npdf << ", for " << nPar
                              << " parameters" << endreq;
            m_PDFCol= new PDF_Data*[m_Npdf];
            for( int pdfs= 0; pdfs<m_Npdf; ++pdfs ) m_PDFCol[pdfs]= 0;
            if( nPar>0 && m_Npdf>0 ) continue;
        } 
        else if ( line == "AXES:" )
        {
            if( !m_PDFAxes )
            {
                log<<MSG::DEBUG<<"Extracting Axes"<<endreq;
                m_PDFAxes= new PDF_Axes(log, dataFile, errflag);
                if( !errflag ) continue;
            }
        } 
        else if ( line.size() > 4 && line.substr(0,4) == "PDF:")
        {
            if( !m_PDFCol ) 
            {
                log << MSG::ERROR << "Line \"#PDFs: <number of PDFs> "
                                  << "Parameters: <number of parameters for PDF function>"
                                  << "\" must come  before 1st PDF";
            } 
            else if (pdf < m_Npdf ) 
            {
                log << MSG::DEBUG << "Extracting PDF Parameters for PDF #"
                                  << pdf << endreq;
                m_PDFCol[pdf++] = new PDF_Data(log, dataFile, m_PDFAxes, nPar, errflag);
                if( !errflag ) continue;
            }
        }
        log << MSG::ERROR << "Incorrect data in file: " << dataFileName << endreq
                          << "Stopping on data line: \"" << line << "\"" << endreq;
        dataFile.close();
        return StatusCode::FAILURE;
    }

    dataFile.close();

    if (pdf < m_Npdf || !m_PDFAxes)
    {
        log << MSG::ERROR << "Missing data in file: " << dataFileName << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}


double CalTkrLikelihoodTool::pdfFCN( double x[2], double parameters[5] ) const
{
    // calculate LogNormal(x[1]*parameters[0]+x[0]):
    // \f[ LogNormal(x, parameters(recEnergy)) \f]
    // if( par[3]==0. ) return 0.; that shouldn't happen
    // \f[ LogNormal(x)= N\,\exp(-\frac{1}{2}(\frac{\ln(1+\tau (x-\mu) \frac{\sinh({\tau}\sqrt{\ln\,4})}{ 2.36\beta{\tau}\,\sqrt{\ln\,4}})}\tau)^2+\tau^2) \f]
    double result = ( x[0] +x[1] * parameters[0] - parameters[2]) / parameters[3];
  
    if( parameters[4]!=0. )
    {
        double val = 1.17741002251547466 * parameters[4]; //sqrt(log(4.))*\tau
        result *= parameters[4] * (sinh(val) / val);
        result  = log(1 + result) / parameters[4]; 
    }
  
    result *= -result * .5;
    result  = parameters[1] * exp(result + parameters[4] * parameters[4]);
    //return isnan(result)?0.:result;
    return result;
}

PDF_Axes::PDF_Axes(MsgStream &log, std::ifstream &dataFile, int &ierr)
{
    m_Sizes      = 0;
    m_BinCenters = 0;
    ierr         = true;

    dataFile >> m_Naxes;

    m_Sizes = new int[m_Naxes];
    m_Bins  = new int[m_Naxes];

    // next are axis sizes
    int ax       = 0;
    int nCenters = 0;
    log << MSG::DEBUG << "Reading Axes ";

    for( ; ax < m_Naxes && !dataFile.eof(); ++ax )
    {
        dataFile >> m_Sizes[ax];
        nCenters += m_Sizes[ax];
        log << MSG::DEBUG << "\t" << m_Sizes[ax];
    }
    
    log << MSG::DEBUG << endreq;
    log << MSG::DEBUG << m_Naxes << " Axes for a total of " << getNbins() << " bins"
                      <<endreq;
    if (ax < m_Naxes) return;
    m_BinCenters = new double[nCenters];

    // next are bin centers
    double* binCenters = m_BinCenters;
  
    ax= 0;
    for( ; ax < m_Naxes && !dataFile.eof(); ++ax )
    {
        int bin = 0;
        log << MSG::DEBUG << "Reading Centers for Axis " << ax << ":" << endreq;
        for( ; bin < m_Sizes[ax] && !dataFile.eof(); ++bin, ++binCenters)
        {
            dataFile >> binCenters[0];
            log << MSG::DEBUG << "(" << bin << ":" << binCenters[0] << ")";
        }
        log << MSG::DEBUG << endreq;
        if( bin != m_Sizes[ax] )
        {
            log << MSG::ERROR << "Wrong Number of values for axis" << ax << endreq;
            return;
        }
    }
    log << MSG::DEBUG << "Read " << (binCenters - m_BinCenters) << " values" << endreq;

    ierr = ax < m_Naxes;

    return;
}

PDF_Axes::~PDF_Axes()
{
    if( m_BinCenters ) delete []m_BinCenters;
    if( m_Sizes )      delete []m_Sizes;
    if( m_Bins )       delete []m_Bins;

    m_Sizes      = 0;
    m_BinCenters = 0;
}

const double* PDF_Axes::getBinCenters(int ax, bool bin) const 
{
    const double *axData= m_BinCenters;
    for( int i= 0; i < ax;  axData += m_Sizes[i++] );
    if(bin) axData += m_Bins[ax];
    return axData;
}

bool PDF_Axes::findBin(int axis, double data) const 
{
    const double* axData = getBinCenters(axis, false);
    int           nabove = m_Sizes[axis]+1;
    int           nbelow = 0;
    int           middle = 0;

    while(nabove - nbelow > 1) 
    {
        middle = (nabove + nbelow) / 2;
        if (data == axData[middle-1]) 
        { 
            nbelow = middle; 
            break; 
        }
        if (data < axData[middle-1])  nabove = middle;
        else nbelow = middle;
    }

    m_Bins[axis] = --nbelow;
    // check if data is not inside axis bounds
    if( nbelow < 0 || nbelow >= m_Sizes[axis] - 1 ) return true;
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
    int id    = 0;
    int width = 1;

    for( int ax = m_Naxes - 1; ax >= 0; --ax )
    {
        id    += width * (m_Bins[ax] + (n & (1 << ax)));
        width *= m_Sizes[ax];
    }
    return id;
}

int PDF_Axes::getNbins(void) const 
{
    int width= 1;
    for( int ax = m_Naxes - 1; ax >= 0; --ax ) width *= m_Sizes[ax];
    return width;
}

PDF_Data::PDF_Data(MsgStream &log, 
                   std::ifstream &dataFile, 
                   PDF_Axes*axes,
                   int npar, 
                   int &ierr):
    m_Npar(npar), m_Axes(axes)
{
    ierr = true;
    if (!m_Axes) return;
    int nBins = m_Axes->getNbins() * m_Npar;
    log << MSG::DEBUG << "Reading " << m_Npar << " parameters per bin for "
                      << nBins/m_Npar << " bins" << endreq;
    m_Data = new float[nBins];
    float* data = m_Data;
    float* end  = m_Data + nBins;
    for( ; data != end && !dataFile.eof(); ++data ) dataFile >> *data;
    log << MSG::DEBUG << "Read " << int(data - m_Data) / m_Npar << " bins" << endreq;
    ierr = data != end;
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

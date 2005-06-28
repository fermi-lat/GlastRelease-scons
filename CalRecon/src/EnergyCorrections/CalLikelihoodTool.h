#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalRecon/ICalEnergyCorr.h"

/**   
* @class CalLikelihoodTool
*
* Performs a binary search for the maximum of a probability density function
* (PDF) returning the most probable value as the event's estimated energy.
* The width of the distribution is used as the error estimate.
* 
* It is used by CalLastLayerLikelihood and CalTkrLikelihood.
*/

class PDF_Axes 
{
  // this class gives the location of phase points at which the PDFs have been
  // estimated. in between those points a linear interpolation is used.
public:
   PDF_Axes(int a, const int *b, const double *c):
     m_BinCenters(c), m_Sizes(b), m_Bins(new int[a]), m_Naxes(a){}
   static PDF_Axes *read(std::ifstream&);
   ~PDF_Axes();

   bool init( const double* );
   int getNbins(void) const;
   int getNaxes(void) const { return m_Naxes; }
   int getSize(int ax) const { return m_Sizes[ax]; }
   const double* getBinCenters(int, bool) const;

   // finds the bin corresponding to the double  on the axis int
   // returns true if not found: point outside PDF phase space
   bool findBin(int, double) const;
   // finds the pahse space bin corresponding to the parameter vector double*
   // returns true if not found: point outside PDF phase space
   bool findBin(double*, int &) const;
   
   // finds the points neighbour to another in the mesh
   int calculateNeigbouringBin(int) const;
private:
    const double *m_BinCenters;
    const int *m_Sizes;
    int *m_Bins;
    
    int m_Naxes;
};

class PDF_Data 
{
  // this class holds the PDF estimates and performs the interpolation
public:
   PDF_Data(float *a, int b, const PDF_Axes *c):
     m_Data(a), m_Axes(c), m_Npar(b){}
   static PDF_Data *read(std::ifstream&, int, const PDF_Axes*);
   ~PDF_Data();

   int getNpar( void ) const { return m_Npar; }
   //preforms the linear interpolation for parameters in const double*
   //returning the estimation in double*
   //returns true if the parameters sit outside the phase space for which the
   //PDFs are defined
   bool evaluateParameters(const double*, double*) const;
private:
    const float *m_Data;
    const PDF_Axes *m_Axes;
    int m_Npar;
};

class CalLikelihoodTool : public AlgTool, virtual public ICalEnergyCorr
{
public:
    //! constructor
    CalLikelihoodTool( const std::string& type, const std::string& name, 
                       const IInterface* parent): AlgTool(type,name,parent){}
    virtual ~CalLikelihoodTool(){}
    
    StatusCode initialize();
    StatusCode finalize();

    //! Energy reconstruction using correlation between energy deposit in the
    // CAL and another observable.
    /* We used the Monte Carlo simulation of the LAT to determine this
    * correlation at discrete incident energies and angles. For each of these,
    * and a parameter \alpha is optimised so as to obtain the best resolution
    * for the distributions:
    *\f[ E_{CAL} + \alpha TkrHits \f]
    * normalised to a probability and with its MPV at the incident energy.
    * All these distributions  can be used to defined a probability density 
    * function. The reconstructed energy for a given event is the one
    * maximising the probability, for a reconstruced direction,
    * CAL raw energy, the TKR vertex, a geometric cut value found by the
    * findGeometricCut, and the correlatad observable.
    *
    *
    * \par The method takes 4 arguments:
    * \param minimum Monte Carlo energy  for which the PDFs are defined
    * \param maximum Monte Carlo energy  for which the PDFs are defined
    * \param CalCluster
    * \param MsgStream implemented in the child class
    *
    * 
    *\return CalCorToolResult with an energy and error estimate.
    *
    *\warning: needs TKR reconstruction
    *\warning: methods setEventPDFdata and setEventPDFparameters should be used
    *\warning:    before calling on calculateEvent
    *
    *\author
    */
           
protected:
    Event::CalCorToolResult *calculateEvent(double, double,
                                            const Event::CalCluster*, 
                                            MsgStream& );
    
    // Sum of the horizontal distance to a crack along the flight path, inside
    // the CAL.
    // This sum is wheighted by the energies inside a layer.
    // It is used to characterise the quality of an event.
    double findGeometricCut( const Point&, const Vector&,
                             const Event::CalCluster* );

    // vertex height is characterised using an int between 0 and 15
    int findTkrVertex( const Point &p ){ return int((p[2]-108.)*3.2e-2); }

    // next 2 methods update the object fields before the use of doEnergyCorr
    void setEventPDFdata( int a ) { m_eventPDF= m_PDFCol[a]; }
    void setEventPDFparameters( double slope, double calE, double xxx ) 
    { m_eventPar[1]= slope; m_eventPar[2]= calE; m_eventPar[3]= xxx; }
    
    // a data set exits which doesn't use any correlation:
    // this function informs the algorithm that sch a set is being used
    bool hasCorrelation(void) const { return m_PDFCol[0]->getNpar()==5; }

    /// input data file containing log normal parameters
    std::string m_dataFile;
    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;

    /// Detector Service
    IGlastDetSvc*     m_detSvc; 

private:
    double& trialEnergy(void) { return m_eventPar[0]; }
    double& calEnergy(void) { return m_eventPar[2]; }
    bool evaluatePDF(double&) const;

    // data tables
    PDF_Axes *m_PDFAxes;
    PDF_Data **m_PDFCol;
    int m_Npdf;

    // CAL origin, tower pitch, and gap between a CAL and its tower,
    // needed for the cuts
    double m_calZorigin;
    double m_towerPitch;
    double m_ratioCDEHeighTowerPitch;
    
    // is the current event PDF, depending on TKR Plane and cut values
    const PDF_Data *m_eventPDF;
    // contains trial energy, slope, CAL Energy and either 0., TKR Hits, or CAL
    // Last Layer Energy
    double m_eventPar[4];
};

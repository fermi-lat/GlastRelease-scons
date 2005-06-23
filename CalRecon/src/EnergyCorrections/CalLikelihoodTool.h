#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "CalRecon/ICalEnergyCorr.h"

class PDF_Axes 
{
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

   bool findBin(int, double) const;
   bool findBin(double*, int &) const;
   int calculateNeigbouringBin(int) const;
private:
    const double *m_BinCenters;
    const int *m_Sizes;
    int *m_Bins;
    
    int m_Naxes;
};

class PDF_Data 
{
public:
   PDF_Data(float *a, int b, const PDF_Axes *c):
     m_Data(a), m_Axes(c), m_Npar(b){}
   static PDF_Data *read(std::ifstream&, int, const PDF_Axes*);
   ~PDF_Data();

   int getNpar( void ) const { return m_Npar; }
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
           
protected:
    double findGeometricCut( const Point&, const Vector&,
                             const Event::CalCluster* );
    int findTkrVertex( const Point &p ){ return int((p[2]-108.)*3.2e-2); }

    void setEventPDFdata( int a ) { m_eventPDF= m_PDFCol[a]; }
    void setEventPDFparameters( double slope, double calE, double xxx ) 
    { m_eventPar[1]= slope; m_eventPar[2]= calE; m_eventPar[3]= xxx; }
    Event::CalCorToolResult *calculateEvent(double, double,
                                            const Event::CalCluster*, 
                                            MsgStream& );
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

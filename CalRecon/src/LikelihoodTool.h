#ifndef __LikelihoodTool_H
#define __LikelihoodTool_H 1

#include "EnergyCorr.h"

/**   
* @class LikelihoodTool
*
* Base class for  energy reconstruction likelihood tools
*
* $Header$
*/

class PDF_Axes {
  private:
    double *m_BinCenters;
    int *m_Sizes;
    int *m_Bins;
    
    int m_Naxes;
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
};

class PDF_Data {
  private:
    float *m_Data;
    int m_Npar;
    const PDF_Axes *m_Axes;
  public:
   PDF_Data(MsgStream&, std::ifstream&, PDF_Axes*, int, int&);
   ~PDF_Data();

   int getNpar( void ) const { return m_Npar; }
   bool evaluateParameters(const double*, double*) const;
};

class LikelihoodTool: public EnergyCorr {
  public:
    typedef enum { cutNOTKRREC= -1, cutMINCALENERGY= -2, cutMAXCALENERGY= -3,
                   cutSLOPE= -4, cutVERTEX= -5, cutPOSITION= -6,
                   cutNOPARAMETERS= -7 }
                 LikelihoodTool_Cuts;
    //! constructor
    LikelihoodTool( const std::string& type, const std::string& name, 
             const IInterface* parent): EnergyCorr(type, name, parent){}
    //! destructor
    ~LikelihoodTool(){}
    StatusCode finalize();
    
    // basic building brick for any PDF
    StatusCode calculateEvent( int, double, double, double*,
                               double*, double&, double& );
    StatusCode readPDFparameters( MsgStream &lm, std::string& );


  private:
    double pdfFCN(double[2], double[5]) const;

    // data tables
    PDF_Axes *m_PDFAxes;
    PDF_Data **m_PDFCol;
    int m_Npdf;
};
#endif

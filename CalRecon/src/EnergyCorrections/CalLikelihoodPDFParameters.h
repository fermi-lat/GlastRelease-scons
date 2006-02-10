#ifndef __CalLikelihoodPDFParameters__H
#define __CalLikelihoodPDFParameters__H
#include <fstream>
class MsgStream;

class PDFGrid
{
  /**
  * @class PDFGrid
  * @author Pol d'Avezac
  * @brief gives the location of phase points at which the PDFs have been
  *    estimated. In other words describes the phase space grid used for
  *    calibration.
  */

  public:
     //!constructor
     PDFGrid(std::ifstream&, MsgStream&);
     //!destructor
     ~PDFGrid();

     // find phase-space point, saved in m_Bins
     bool initialise( const double[]);

     // returns the number of poins on the grid
     int getNbins(void) const;
     // returns the number of axes for the grid
     int getNaxes(void) const { return m_Naxes; }
     // returns the number of points on axis ax
     int getSize(int ax) const { return m_Sizes[ax]; }

     double minTrialEnergy(void) const { return m_BinCenters[0]; }
     double maxTrialEnergy(void) const { return m_BinCenters[m_Sizes[0]-1]; }
    
     // finds the bin corresponding to the double  on the axis int
     // returns true if not found: point outside PDF phase space
     bool findBin(int, double) const;
     // finds the pahse space bin corresponding to the parameter vector double*
     // returns true if not found: point outside PDF phase space
     bool findBin(double*, int &) const;
     // finds the points neighbour to another in the mesh
     int calculateNeigbouringBin(int) const;

     double getBinCenter(int ax, int bin) const;
     const double* getBinCenters(int, bool) const;

     // returns wether the object has been initialised
     bool operator!(void){ return !m_BinCenters; }
  private:
     // contains bin grid: (energy, angle) values at which the PDFParameterss were
     // calculated 
     const double *m_BinCenters;
     // contains grid widths
     const int *m_Sizes;

     // used for interpolation: it defines one pont inside the grid
     int *m_Bins;
     // number of axes in the grid
     int m_Naxes;
};

class PDFObject
{
  /**
  * @class PDFObject
  * @author Pol d'Avezac
  * @brief wrapper class used by PDFAnswer
  */
  public:
   PDFObject(){}
   virtual ~PDFObject(){}
   
   // returns wether the object was initialised
   virtual bool operator!(void) const { return true; }
};

class PDFParameters: public PDFObject
{
  /**
  * @class PDFParameters
  * @author Pol d'Avezac
  * @brief contains PDFFunction parameters for each point on a (energy, angle to
  *    the vertical) grid described by the PDFGrid object.
  *    It performs a linear interpolation to estimate the parameters between
  *    points on the grid.
  */
  public:
     //!constructor
     PDFParameters(std::ifstream&, PDFGrid*, int nPar, MsgStream &log);
     //!destructor
     ~PDFParameters();

     bool initialise( const double par[] ){ return m_Grid->initialise(par); }
     bool interpolation(const double[], double[], bool);

     // number of parameters at each grid point
     int getNpar( void ) const { return m_Npar; }
     const PDFGrid* getGrid(void) const { return m_Grid; }

     double minTrialEnergy(void) const { return m_Grid->minTrialEnergy(); }
     double maxTrialEnergy(void) const { return m_Grid->maxTrialEnergy(); }

     // returns wether the object has been initialised
     bool operator!(void) const { return !m_Data; }
  private:
     // number of parameters at each grid point
     int m_Npar;
     // parameters at each grid point
     const double *m_Data;

     // describes the grid
     PDFGrid *m_Grid;
};
#endif

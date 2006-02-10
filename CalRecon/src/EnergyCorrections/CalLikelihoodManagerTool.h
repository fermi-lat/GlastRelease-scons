#ifndef __CalLikelihoodManagerTool__H
#define __CalLikelihoodManagerTool__H
#include "GaudiKernel/AlgTool.h"
#include "CalRecon/ICalEnergyCorr.h"
#include "CalLikelihoodPDFParameters.h"

#include "Event/Digi/TkrDigi.h"

#include <fstream>
#include <vector>
/**   
* @class  CalLikelihoodManagerTool
* @author Pol d'Avezac
*
* @brief Tool to correct energy using maximum likelihood
*
* It uses a number of different classes:
*   PDFFunction: describe either the form of the likelihood or of cuts
*   PDFParameters: contains data necessary for PDFFunction. This data was calculated
*     for fixed values of energy and angle to the vertical. This object provides
*     an interpolation method.
*   PDFGrid: contains the (energy, angle to the vertical) grid used by PDFParameters
*
*   PDFAnswer: base class for all those providing an energy estimation.
*     It provides the mecanisms for choosing the best estimation and accessing
*     it.
*   PDFVertexArray: contains all PDFParameters for 1 cut value, given 1 cut definition
*   PDFCutArray: contains all PDFVertexArrays for every cuts in 1 cut definition
*     It contains the PDFVertexArray and PDFFunction for defining this cut and
*     provides a method returning wheter or not a PDFVertexArray estimation in within
*     its cut boundaries.
*   CalLikelihoodManagerTool: contains all PDFCutArrays. There are at the
*     moment 3 such tools for 2 different cut definitions and 3 (energy, angle)
*     range. It provides access to the Gaudi environment for all other classes.
*
* $Header:
*/

typedef std::vector<PDFObject*> PDFVect;
typedef std::vector<PDFObject*>::iterator PDFVectItr;
typedef std::vector<PDFObject*>::const_iterator PDFVectConstItr;

class PDFAnswer: public PDFObject, public PDFVect
{
  /**
  *  @class PDFAnswer
  *  @author Pol d'Avezac
  *  @brief base class for all classes retruning an energy estimate.
  
  *  It implements the recursive calculations on all objects in PDFVect and
  *  holds the one with the best energy estimate in its m_Answer field.
  *  It uses that to implement to recursive access to the best energy estimate
  */

  public:
    typedef enum { kNotCalculated= int(0),
                   kUnphysical= int(0),
                   kNoCuts= int(2),
                   kOutsideCut= int(3),
                   kInsideCut= int(4) } Status_t;
    //!constructor
    PDFAnswer(): PDFObject(), PDFVect(0), m_Answer(0){}
    //!destructor
    ~PDFAnswer(){}

    virtual void clear(void);

    // returns an estimation of the energy in the range given by the first 2
    // doubles
    virtual bool evalEnergy(double, double, const Event::CalCluster*,
                                            const Event::TkrVertex*);
    // returns an estimation of the error
    virtual bool evalError(Event::TkrVertex *p) { return m_Answer?m_Answer->evalError(p): true; }
    
    // returns vertex altitude
    int getTkrPlane(const Event::TkrVertex *vertex) const
    { return int(3.2e-2*(vertex->getPosition().z()-108.)); }

    // returns estimated energy
    virtual double getEnergy(void) const { return m_Answer?m_Answer->getEnergy(): 0; }
    // returns estimated probability
    virtual double getProbability(void) const { return m_Answer?m_Answer->getProbability(): 0; }
    // returns estimated errror
    virtual double getError(void) const { return m_Answer?m_Answer->getError(): 0; }
    // returns estimated lower or upper range error
    virtual double getError(bool a) const { return m_Answer?m_Answer->getError(a): 0; }

    // returns object within PDFVect with the best energy estimate
    virtual const PDFAnswer *getAnswer(bool a= false) const
    { if( a ) return m_Answer?m_Answer->getAnswer(true): 0; else return m_Answer; }
    void setAnswer(PDFAnswer *a) { m_Answer= a; }
    // returns status: inside or ouside cuts boundaries
    virtual int getStatus() const { return m_Answer?m_Answer->getStatus(): kNotCalculated; }

    int getAnswerID(void) const;
    // used to find best energy estimate
    bool operator<(const PDFAnswer &) const;
    bool operator!(void) const { return PDFVect::empty(); }
  protected:
    // used for initialisation from file
    static bool checkLine(std::string &);
    static bool checkField(std::string &, const char[]);
  private:
    PDFAnswer *m_Answer;
};

class PDFFunction;
class PDFCutArray;
class CalLikelihoodManagerTool;
class IGlastDetSvc;
class IDataProviderSvc;

class PDFVertexArray: public PDFAnswer
{
  /**
  *  @class PDFVertexArray
  *  @author Pol d'Avezac
  *  @brief class containing likelihoods for 1 cut value and all TkrVertex
  *  altitudes
  *
  *  This is the first object in the PDFAnswer recursion
  */
  public:
      //! constructor
      PDFVertexArray(std::ifstream&, PDFCutArray*, PDFFunction*, MsgStream&);
    //!destructor
      ~PDFVertexArray(){}
      
      void clear(void) { PDFAnswer::clear(); m_Fcn= 0; m_Tool= 0; }

      // returns an estimation of the energy in the range given by the first 2
      // doubles
      bool evalEnergy(double, double, const Event::CalCluster*,
                                      const Event::TkrVertex*);
      // returns an estimation of the error
      bool evalError(Event::TkrVertex *p);

      // returns estimated energy
      double getEnergy(void) const { return m_MPV[0]; }
      // returns estimated probability
      double getProbability(void) const { return m_MPV[1]; }
      // returns estimated error
      double getError(void) const { return (m_FWHM[1]+m_FWHM[0])*.5; }
      // returns estimated lower or upper range error
      double getError(bool a) const { return m_FWHM[a]; }

      const PDFAnswer *getAnswer(bool= false) const { return this; }
      int getStatus() const { return m_status; }
  private:
      // defines the likelihood form.
      PDFFunction *m_Fcn;

      // pointer to its container
      PDFCutArray *m_Tool;

      // best energy and its probability
      double m_MPV[2];
      // lower and upper range error
      double m_FWHM[2];
      // status: inside or outside cut boundaries
      Status_t m_status;
};

class PDFCutArray: public PDFAnswer
{
  /**   
  * @class PDFCutArray
  * @author Pol d'Avezac
  *
  * @brief Tool returning the best energy estimation from all cut values for 1
  *     cut definition
  *
  * The cut definition uses : PDFVertexArray m_Cuts
  *                           PDFFunction m_CutsFcn
  * The likelihoods with parameters in the vector of PDFAnswers 
  * The tool works in the range described by the PDFGrid range
  * PDFCutArray holds 1 Likelihood per set of cut and altitude.
  * Cut sets are described using PDFParameters and PDFFunction classes as well.
  * Likelihoods for 1 cut value and all Tkr vertexes are regrouped in 
  * PDFVertexArray
  * 
  *
  * $Header:
  */
  public:
      //! constructor
      PDFCutArray(std::string, const CalLikelihoodManagerTool*, MsgStream&);
      //! destructor
      ~PDFCutArray(){}

      void clear(void);

      // returns an estimation of the energy in the range given by the first 2
      // doubles
      bool evalEnergy(double, double, const Event::CalCluster*,
                                      const Event::TkrVertex*);
      // returns status of PDFVertexArray for vertex altitude int
      Status_t evalStatus(int, PDFVertexArray*);

      // number of different cut ranges
      int getNarray(void){ return m_Narray; }
      // number of PDFParameters
      int getNpdf(void){ return m_Npdf; }
      PDFGrid* getGrid(void) { return m_Grid; }

      PDFVertexArray* getCuts(void) { return m_Cuts; }
      PDFFunction* getLikelihoodFcn(void) { return m_LikelihoodFcn; }
      PDFFunction* getCutsFcn(void) { return m_CutsFcn; }
  private:
      // number of different cut ranges
      int m_Narray;
      // number of PDFParameters
      int m_Npdf;

      // Defines the (energy, angle) grid at which PDFParameters values were
      // created used for interpolating of these values
      PDFGrid *m_Grid;

      // Defines the fucntion describing the form of the likelihood
      PDFFunction *m_LikelihoodFcn;
      // Defines the fucntion describing the form of the cuts
      PDFFunction *m_CutsFcn;
      // Contains PDFParameters for all vertex altitudes
      PDFVertexArray *m_Cuts;
};

class CalLikelihoodManagerTool: public AlgTool, virtual public ICalEnergyCorr,
                                public PDFAnswer{
  public:
    //! constructor
    CalLikelihoodManagerTool(const std::string&, const std::string&, 
                             const IInterface*);
    ~CalLikelihoodManagerTool(){}
    
    StatusCode initialize();
    StatusCode finalize();

    //! Energy reconstruction using correlation between energy deposit in the
    // CAL and another observable.
    /* We used the Monte Carlo simulation of the LAT to determine this
    * correlation at discrete incident energies and angles. For each of these,
    * and a parameter \alpha is optimised so as to obtain the best resolution
    * for the distributions:
    *\f[ CalEnergyRaw + \alpha CAlELayer7 + \beta TkrSumHits \f]
    * normalised to a probability and with its MPV at the incident energy.
    * All these distributions  can be used to defined a probability density 
    * function. The reconstructed energy for a given event is the one
    * maximising the probability, for a reconstruced direction,
    * CalEnergyRaw, ...
    *
    *
    * \par The method takes 2 arguments:
    * \param Event::CalCluster*
    * \param Event::TkrCluster*
    *
    *\return CalCorToolResult with an energy and error estimate.
    *
    *\warning: needs TKR reconstruction
    *
    *\author Pol d'Avezac
    */
    Event::CalCorToolResult* doEnergyCorr(Event::CalClusterCol*, 
                                          Event::TkrVertex*);

    bool push_back(const char[], MsgStream&);
    void getParameters(std::map<double*,std::string>&, MsgStream&) const;
    const Event::TkrDigiCol *getTkrDigiCol(void) const;
    
  private:
    double m_minTrialEnergy;
    double m_maxTrialEnergy;
    double m_minTkr1ZDir;
           

    //! input data files containing log normal parameters
    StringArrayProperty   m_parameterFiles ;

    //! Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;

    //! Detector Service
    IGlastDetSvc*     m_detSvc; 
};
#endif

#ifndef __CalLikelihoodPDFFunctions__H
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#define __CalLikelihoodPDFFunctions__H

namespace Event {
  class CalCluster;
  class TkrVertex;
}
class PDFParameters;
class CalLikelihoodManagerTool;
class MsgStream;

class PDFFunction
{
   /**
   * @class PDFFunction
   * @author Pol d'Avezac
   * @brief base class for describing the shape of either likelihood or cuts.
   *
   *  This class uses parameters defined in PDFParameters objects to return a
   *  cut value or a probability. The actual shape definition is done in derived
   *  classes. It provides with methods:
   *    - value: returning the function value at a given point
   *    - getMPV: returning the MPV and its value in a given MC Energy range
   *    - getFWHM: returning the FWHM for a previously calculated MPV.
   */
   
   public:
     //!constructor
     PDFFunction(int nIntPar, int nFcnPar, int nFcnArg):
       m_Nstep(15),
       m_NintPar(nIntPar), m_IntPar(nIntPar?new double[nIntPar]:0),
       m_NfcnPar(nFcnPar), m_FcnPar(nFcnPar?new double[nFcnPar]:0),
       m_NfcnArg(nFcnArg), m_FcnArg(nFcnArg?new double[nFcnArg]:0){}
     //!destructor
     virtual ~PDFFunction(){ clear(); }

     //implemented and used in derived PDF...Cuts functions.
     virtual bool unPhysical(double){ return false; }
     //implemented in derived classes only
     virtual bool value(double ans[1]) const= 0;
     virtual void setEvt(const Event::CalCluster *cal, 
                         const Event::TkrVertex*)= 0;

     bool eval(PDFParameters*, double[1], bool= true);

     bool getMPV(PDFParameters*, double, double, double[2]);
     bool getFWHM(PDFParameters*, const double[2], double[2]);

     void setNfunctionParameters(int);
     void setNfunctionArgumentents(int);
     void setNinterpolationParameters(int);
     void clear(void);

     int getNfunctionParameters(void) const { return m_NfcnPar; }
     int getNfunctionArguments(void) const { return m_NfcnArg; }
     int getNinterpolationParameters(void) const { return m_NintPar; }
     const double* getInterpolationParameters(void) const { return m_IntPar; }
     const double* getFunctionParameters(void) const { return m_FcnPar; }
     const double* getFunctionArguments(void) const { return m_FcnArg; }
     double* getInterpolationParameters(void) { return m_IntPar; }
     double* getFunctionArguments(void) { return m_FcnArg; }

     // The MPV is looked for using a divide-and-conquer method
     // This is the number of evaluations for each recursion
     int getNstep( void ) const { return m_Nstep; }
     void setNstep(int a) { m_Nstep= a; }

     // getMPV current estimation of the energy
     double& trialEnergy(void) { return m_IntPar[0]; }
   private:
     // CalEnergyRaw value, used in getFWHM
     double& firstArg(void) { return m_FcnArg[0]; }
     // interval width for the divide-and-conquer getMPV method
     double delta(double a, double b) const { return (b-a)/m_Nstep; }

     int m_Nstep;

     int m_NintPar;
     double *m_IntPar;

     int m_NfcnPar;
     double *m_FcnPar;

     int m_NfcnArg;
     double *m_FcnArg;
};

class PDFLikelihood: public PDFFunction {
  /**
   * @class PDFLikelihood
   * @author Pol d'Avezac
   * @brief class describing the Likelihood shape
   */
 public:
  PDFLikelihood(const CalLikelihoodManagerTool *man, MsgStream&);
  void setEvt(const Event::CalCluster*, const Event::TkrVertex*);
  bool value(double[1]) const;
  
  bool unPhysiscal(double energy) const { return calEnergyRaw()>energy; }
  
  // next functions are really for comfort
  double tkr1ZDir(void) const { return getInterpolationParameters()[1]; }
  double& tkr1ZDir(void) { return getInterpolationParameters()[1]; }
  double calEnergyRaw(void) const { return getFunctionArguments()[0]; }
  double& calEnergyRaw(void) { return getFunctionArguments()[0]; }
  double calELayer7(void) const { return getFunctionArguments()[1]; }
  double& calELayer7(void) { return getFunctionArguments()[1]; }
  double tkrSumHits(void) const { return getFunctionArguments()[2]; }
  double& tkrSumHits(void) { return getFunctionArguments()[2]; }
  
  double occupancy(void) const { return getFunctionParameters()[0]; }
  double calELayer7Alpha(void) const { return getFunctionParameters()[1]; }
  double tkrSumHitsAlpha(void) const { return getFunctionParameters()[2]; }
  double norm(void) const { return getFunctionParameters()[3]; }
  double mpv(void) const { return getFunctionParameters()[4]; }
  double sigma(void) const { return getFunctionParameters()[5]; }
  double tau(void) const { return getFunctionParameters()[6]; }
 private:
  const CalLikelihoodManagerTool *m_manager;
};

class PDFLowEnergyCuts: public PDFFunction {
  /**
   * @class PDFLikelihood
   * @author Pol d'Avezac
   * @brief class describing the cut boundaries for:
   *    -Low energy photons
   *    -Hig henergy, high incidence photons
   */
 public:
  typedef enum { cZECNTR_MIN= int(0), 
		 cZECNTR_MAX= int(1) } Cuts_t;
  PDFLowEnergyCuts(const CalLikelihoodManagerTool*, MsgStream&);
  
  void setEvt(const Event::CalCluster*, const Event::TkrVertex*);
  bool value(double[1]) const;
  
  // next functions are really for comfort
  double tkr1ZDir(void) const { return getInterpolationParameters()[1]; }
  double& tkr1ZDir(void) { return getInterpolationParameters()[1]; }
  double calZEcntr(void) const { return getFunctionArguments()[0]; }
  double& calZEcntr(void) { return getFunctionArguments()[0]; }
  double geometricCut(void) const { return getFunctionArguments()[1]; }
  double& geometricCut(void) { return getFunctionArguments()[1]; }
  
  // calculates the geometricCut parmeter's value
  double geometricCut(const Event::CalCluster*,
		      const Event::TkrVertex*) const;
 private:
  const CalLikelihoodManagerTool *m_manager;
  double m_calZorigin;
  double m_towerPitch;
  double m_ratioCDEHeighTowerPitch;
};

class PDFHighEnergyCuts: public PDFFunction {
  /**
   * @class PDFLikelihood
   * @author Pol d'Avezac
   * @brief class describing the cut boundaries for:
   *    -Hig henergy, low incidence photons
   */
 public:
  PDFHighEnergyCuts(const CalLikelihoodManagerTool*, MsgStream&);
  
  void setEvt(const Event::CalCluster*, const Event::TkrVertex*);
  bool value(double[1]) const;
  
  // next functions are really for comfort
  double tkr1ZDir(void) const { return getInterpolationParameters()[1]; }
  double& tkr1ZDir(void) { return getInterpolationParameters()[1]; }
  double calZEcntr(void) const { return getFunctionArguments()[0]; }
  double& calZEcntr(void) { return getFunctionArguments()[0]; }
  double calTwrEdgeCntr(void) const { return getFunctionArguments()[1]; }
  double& calTwrEdgeCntr(void) { return getFunctionArguments()[1]; }
  double calELayer7(void) const { return getFunctionArguments()[2]; }
  double& calELayer7(void) { return getFunctionArguments()[2]; }
  
  // calculates the CalTwrEdgeCntr parmeter's value
  // similar to the Merit's synonym
  double calTwrEdgeCntr(const Event::CalCluster *cluster)const;
 private:
  const CalLikelihoodManagerTool *m_manager;
  typedef enum { cZECNTR_TOP      = int(0), cZECNTR_CORE     = int(1),
		 cZECNTR_BOTTOM   = int(2), cE7          = int(3),
		 cTOWER_EDGE  = int(4), cTOWER_BORDER= int(5),
		 cTOWER_CENTER= int(6) } Cuts_t;
  double m_towerPitch;
};
#endif

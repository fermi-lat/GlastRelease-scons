/**
 * @file Fitter.h
 * @brief Base class for fitting functions to Root distribtions.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_Fitter_h
#define RootAnalysis_irfAnalysis_Fitter_h

class TH1F;

/**
 * @class Fitter
 * @brief Fitter base class.
 * @author J. Chiang
 *
 * $Header$
 */

class Fitter {

public:

   virtual ~Fitter(){}

   virtual void applyFit(TH1F *) = 0;

   virtual void setBounds(const std::vector<double> &,
                          const std::vector<double> &) = 0;

protected:

   Fitter(){}

};

#endif // RootAnalysis_irfAnalysis_Fitter_h

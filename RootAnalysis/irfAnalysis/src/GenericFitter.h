/**
 * @file GenericFitter.h
 * @brief Fit a linear model, y = [0]*x + [1]
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_GenericFitter_h
#define RootAnalysis_irfAnalysis_GenericFitter_h

#include "Fitter.h"

/**
 * @class GenericFitter
 * @brief Derived class of Fitter.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class GenericFitter : public Fitter {

public:

   GenericFitter(const std::string &usrFunction,
                 const std::vector<double> &params, 
                 const std::string &outputFile = "", 
                 int maxTrys=5);

   virtual ~GenericFitter();

   virtual void applyFit(TH1 *);

private:

   int m_maxTrys;

};

#endif // RootAnalysis_irfAnalysis_GenericFitter_h

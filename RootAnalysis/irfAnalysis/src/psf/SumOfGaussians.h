/**
 * @file SumOfGaussians.h
 * @brief Fit a sum of two Gaussians to a TH1F object.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_SumOfGaussians_h
#define RootAnalysis_irfAnalysis_SumOfGaussians_h

#include "Fitter.h"

/**
 * @class SumOfGaussians
 * @brief Derived class of Fitter.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class SumOfGaussians : public Fitter {

public:

   SumOfGaussians(const std::string &outputFile = "", int maxTrys=5);

   virtual ~SumOfGaussians();

   virtual void applyFit(TH1 *);

private:

   int m_maxTrys;

};

#endif // RootAnalysis_irfAnalysis_SumOfGaussians_h

/**
 * @file SumOfGaussians.h
 * @brief Fit a sum of two Gaussians to a TH1F object.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_SumOfGaussians_h
#define RootAnalysis_irfAnalysis_SumOfGaussians_h

#include <vector>

#include "Fitter.h"

class TF1;

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

   SumOfGaussians(int maxTrys=5);

   virtual ~SumOfGaussians();

   virtual void applyFit(TH1F *);

   void setBounds(const std::vector<double> &lower,
                  const std::vector<double> &upper);

//    void setParameters(const std::vector<double> &params);
//    void getParameters(std::vector<double> &params) const;

private:

   TF1 * m_func;
   int m_maxTrys;
};

#endif // RootAnalysis_irfAnalysis_SumOfGaussians_h

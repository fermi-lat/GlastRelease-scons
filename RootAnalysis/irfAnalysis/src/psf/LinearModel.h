/**
 * @file LinearModel.h
 * @brief Fit a linear model, y = [0]*x + [1]
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_LinearModel_h
#define RootAnalysis_irfAnalysis_LinearModel_h

#include "Fitter.h"

/**
 * @class LinearModel
 * @brief Derived class of Fitter.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class LinearModel : public Fitter {

public:

   LinearModel(const std::string &outputFile = "", int maxTrys=5);

   virtual ~LinearModel();

   virtual void applyFit(TH1 *);

private:

   int m_maxTrys;

};

#endif // RootAnalysis_irfAnalysis_LinearModel_h

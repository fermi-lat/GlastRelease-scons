/**
 * @file PsfModel.h
 * @brief Fit the ad hoc model, 
 *    [0]*x*exp(-0.5*pow((x - [1])/[2], 2))*exp(-sqrt(x/[3]))
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_PsfModel_h
#define RootAnalysis_irfAnalysis_PsfModel_h

#include "Fitter.h"

/**
 * @class PsfModel
 * @brief Derived class of Fitter.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class PsfModel : public Fitter {

public:

   PsfModel(const std::string &outputFile = "", int maxTrys=5);

   virtual ~PsfModel();

   virtual void applyFit(TH1 *);

private:

   int m_maxTrys;

};

#endif // RootAnalysis_irfAnalysis_PsfModel_h

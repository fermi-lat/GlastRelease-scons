/**
 * @file Fitter.h
 * @brief Base class for fitting functions to Root distribtions.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_Fitter_h
#define RootAnalysis_irfAnalysis_Fitter_h

#include <vector>
#include <string>

class TH1;
class TFile;
class TTree;
class TF1;

/**
 * @class Fitter
 * @brief Fitter base class.
 * @author J. Chiang
 *
 * $Header$
 */

class Fitter {

public:

   Fitter(const std::string &outputFile);

   virtual ~Fitter();

   virtual void applyFit(TH1 *) = 0;

   virtual void setBounds(const std::vector<double> &,
                          const std::vector<double> &);

   virtual void setParameters(const std::vector<double> &);

   virtual void getParameters(std::vector<double> &) const;

   virtual void writeFitPars();

protected:

   TFile * m_file;
   TTree * m_tree;
   TF1 * m_func;
   std::vector<Double_t> m_params;

};

#endif // RootAnalysis_irfAnalysis_Fitter_h

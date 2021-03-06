/**
 * @file MakeDists.h
 * @brief Interface for creating distributions.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_MakeDists_h
#define RootAnalysis_irfAnalysis_MakeDists_h

#include <vector>

#include "PSF.h"

class TF1;
class Fitter;

/**
 * @class MakeDists
 * @brief Create distributions using Root's project method.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class MakeDists : public PSF {

public:

   MakeDists(const std::string &summaryFile, bool makeProfile=false,
             bool useLogScale=false);

//   ~MakeDists(){}
   
   void project(const std::string &branchName, 
                double xmin, double xmax, int nbins=50,
                Fitter *fitter=0);

   void draw(const std::string &ps_filename, bool logy=false, Fitter* fit=0);

   void setEnergyScaling(std::string scalingFunction,
                         const std::vector<double> &params,
                         bool isThin=true);

   void addEnergyScaling(const std::string &rootFile, 
                         const std::string &treeName);

   void addCutInfo(const std::string &rootFile,
                   const std::string &treeName);

private:

   bool m_makeProfile;
   bool m_useLogScale;

   int m_nbins;

   std::string m_branchName;

   std::string m_scalingFunction;
   std::vector<double> m_params;
   bool m_isThin;

   TF1 * m_energyScale;

   void addThetaGrid(const std::string &rootFile);
   void addEnergyGrid(const std::string &rootFile);

   void applyEnergyScaling();
   void replaceVariables();
   void modifyBranchName();
   void replace_substring(std::string & expression, const std::string &oldstr,
                          const std::string &newstr);

};

#endif // RootAnalysis_irfAnalysis_MakeDists_h

/**
 * @file MakeDists.h
 * @brief Interface for creating distributions.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef RootAnalysis_irfAnalysis_MakeDists_h
#define RootAnalysis_irfAnalysis_MakeDists_h

#include "PSF.h"

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

//    MakeDists(const std::string &summaryFile, bool applyEnergyScaling=false,
//              bool makeProfile=false);
   MakeDists(const std::string &summaryFile,
             bool makeProfile=false);

//   ~MakeDists(){}
   
   void project(const std::string &branchName, 
                double xmin, double xmax, int nbins=50,
                Fitter *fitter=0);


   void draw(const std::string &ps_filename,  bool logy=false, Fitter* fit=0);

//    void setEnergyScaling(const std::string &scalingFunction);

//    void addEnergyScaling(const std::string &rootFile, 
//                          const std::string &treeName);

   void addCutInfo(const std::string &rootFile,
                   const std::string &treeName);

private:

   bool m_makeProfile;

   int m_nbins;

//    double energyScaling(double energy);

};

#endif // RootAnalysis_irfAnalysis_MakeDists_h

/**
 * @file SumOfGaussians.cxx
 * @brief Implementation for fitting sum of two Gaussians to a Root
 * histogram.
 * @author J. Chiang
 *
 * $Header$
 */

#include <sstream>

#include "TF1.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"

#include "SumOfGaussians.h"

SumOfGaussians::SumOfGaussians(const std::string &outputFile, int maxTrys) 
   : Fitter(outputFile), m_maxTrys(maxTrys) {
// Use a more natural formulation of the Gaussian function than Root's
// built-in "gaus" function.  The ordering for a single Gaussian is 
// [0] = integral, [1] = x0, [2] = sigma.
   std::ostringstream gaussSum;
   gaussSum << "[0]/sqrt(2.)/[2]*exp(-0.5*pow((x - [1])/[2], 2.)) "
            << "+ [3]/sqrt(2.)/[5]*exp(-0.5*pow((x - [4])/[5], 2.))";

   m_func = new TF1("gaussSum", gaussSum.str().c_str());

   if (m_tree) {
// Set up the branches in the output tree.
      m_params.resize(6);
      m_tree->Branch("norm1" , &m_params[0], "norm1/D");
      m_tree->Branch("x01"   , &m_params[1], "x01/D");
      m_tree->Branch("sigma1", &m_params[2], "sigma1/D");
      m_tree->Branch("norm2" , &m_params[3], "norm2/D");
      m_tree->Branch("x02"   , &m_params[4], "x02/D");
      m_tree->Branch("sigma2", &m_params[5], "sigma2/D");
   }
}

SumOfGaussians::~SumOfGaussians() {
}

void SumOfGaussians::applyFit(TH1 * h) {
// Unfortunately, this initialization must be fairly context-specific
// to work well, here fitting the log10(Tkr1[Theta,Phi]Err)
// distributions.
   Double_t scale = h->Integral();
   Double_t mean = h->GetMean();
   Double_t rms = h->GetRMS();

   std::vector<Double_t> params;
   params.push_back(scale);
   params.push_back(mean);
   params.push_back(rms);
   params.push_back(scale/2.);
   params.push_back(-1.8);
   params.push_back(rms);

   m_func->SetParameters(&params[0]);

   int fitTrys = 0;
// Do the fit in quiet mode.
   while (h->Fit("gaussSum", "Q") && fitTrys++ < m_maxTrys);
}

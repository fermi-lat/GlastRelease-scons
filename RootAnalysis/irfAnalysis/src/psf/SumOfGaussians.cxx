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

#include "SumOfGaussians.h"
#include <cassert>

SumOfGaussians::SumOfGaussians(int maxTrys) : m_func(0), m_maxTrys(maxTrys) {
// Use a more natural formulation of the Gaussian function than Root's
// built-in "gaus" function.  The ordering for a single Gaussian is 
// [0] = integral, [1] = x0, [2] = sigma.
   std::ostringstream gaussSum;
   gaussSum << "[0]/sqrt(2.)/[2]*exp(-0.5*pow((x - [1])/[2], 2.)) "
            << "+ [3]/sqrt(2.)/[5]*exp(-0.5*pow((x - [4])/[5], 2.))";

   m_func = new TF1("gaussSum", gaussSum.str().c_str());
}

SumOfGaussians::~SumOfGaussians() {
   delete m_func;
}

void SumOfGaussians::applyFit(TH1F * h) {
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
   while (h->Fit("gaussSum") && fitTrys++ < m_maxTrys);
}

void SumOfGaussians::setBounds(const std::vector<double> &lower,
                               const std::vector<double> &upper) {
   assert(lower.size() == upper.size());
   for (unsigned int i = 0; i < lower.size(); i++) {
      m_func->SetParLimits(i, lower[i], upper[i]);
   }
}

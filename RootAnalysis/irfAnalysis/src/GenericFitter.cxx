/**
 * @file GenericFitter.cxx
 * @brief Provides a convenient interface to Root's Fit method.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cassert>
#include <algorithm>
#include <sstream>

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"

#include "GenericFitter.h"

GenericFitter::GenericFitter(const std::string &userFunction,
                             const std::vector<double> &params,
                             const std::string &outputFile, 
                             int maxTrys) 
   : Fitter(outputFile), m_maxTrys(maxTrys) {

   m_func = new TF1("userFunction", userFunction.c_str());

   unsigned int npar = m_func->GetNpar();
   m_params.resize(npar);
   assert(npar >= params.size());
   std::copy(params.begin(), params.end(), m_params.begin());

// Initialize parameters; subsequent fits will start with parameters
// from the previous fit.
   m_func->SetParameters(&m_params[0]);

   if (m_tree) {
// Set up the branches in the output tree.
      for (unsigned int i = 0; i < npar; i++) {
         std::ostringstream parName;
         parName << "p" << i;
         m_tree->Branch(parName.str().c_str(), &m_params[i], 
                        (parName.str() + "/D").c_str());
      }
   }
}

GenericFitter::~GenericFitter() {
}

void GenericFitter::applyFit(TH1 * h) {
   int fitTrys = 0;

// Do the fit in quiet mode.
   while (h->Fit("userFunction", "Q") && fitTrys++ < m_maxTrys);
}

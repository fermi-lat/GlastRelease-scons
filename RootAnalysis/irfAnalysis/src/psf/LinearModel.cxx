/**
 * @file LinearModel.cxx
 * @brief Implementation for fitting y = [0]*x +[1] to a Root profile.
 * @author J. Chiang
 *
 * $Header$
 */

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"

#include "LinearModel.h"

LinearModel::LinearModel(const std::string &outputFile, int maxTrys) 
   : Fitter(outputFile), m_maxTrys(maxTrys) {

   m_func = new TF1("linearModel", "[0]*x + [1]");

   if (m_tree) {
// Set up the branches in the output tree.
      m_params.resize(2);
      m_tree->Branch("slope" , &m_params[0], "slope/D");
      m_tree->Branch("intercept" , &m_params[1], "intercept/D");
   }
}

LinearModel::~LinearModel() {
}

void LinearModel::applyFit(TH1 * h) {
   int fitTrys = 0;
// Do the fit in quiet mode.
   while (h->Fit("linearModel", "Q") && fitTrys++ < m_maxTrys);
}

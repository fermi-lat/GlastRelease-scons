/**
 * @file Fitter.cxx
 * @brief Implementation for Fitter base class.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cassert>
#include <algorithm>

#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

#include "Fitter.h"

Fitter::Fitter(const std::string &outputFile) 
   : m_file(0), m_tree(0), m_func(0) {

   if (outputFile != "") {
// Write to the Root output directory.
      std::string outFile = ::getenv("output_file_root")
         + std::string("/") + outputFile;

// Create the new Root file.
      m_file = new TFile(outFile.c_str(), "recreate");

// Create the TTree.
      m_tree = new TTree("fitParams", "Fit parameters");

// Sub-class constructor will create the branches.
   }
}

Fitter::~Fitter() {
   if (m_tree) {
      m_file->cd();
      m_tree->Write();
      m_tree->Print();
   }
//   delete m_func;
   delete m_tree;
   delete m_file;
}

void Fitter::setBounds(const std::vector<double> &lower,
                       const std::vector<double> &upper) {
   assert(lower.size() == upper.size());
   assert(lower.size() <= static_cast<unsigned int>(m_func->GetNpar()));
   for (unsigned int i = 0; i < lower.size(); i++) {
      m_func->SetParLimits(i, lower[i], upper[i]);
   }
}

void Fitter::setParameters(const std::vector<double> &params) {
   assert(static_cast<unsigned int>(m_func->GetNpar()) == params.size());
// This copy may not be necessary.
   std::vector<Double_t> my_params(params.size());
   std::copy(params.begin(), params.end(), my_params.begin());
   m_func->SetParameters(&my_params[0]);
}

void Fitter::getParameters(std::vector<double> &params) const {
   std::vector<Double_t> my_params(m_func->GetNpar());
   m_func->GetParameters(&my_params[0]);
   params.resize(m_func->GetNpar());
   std::copy(my_params.begin(), my_params.end(), params.begin());
}

void Fitter::writeFitPars() {
   if (m_tree) {
      assert(m_params.size() == static_cast<unsigned int>(m_func->GetNpar()));
      m_func->GetParameters(&m_params[0]);
      m_tree->Fill();
   }
}

double Fitter::integral(double a, double b, double eps) {
// Use TF1's numerical integrator (dgauss).  The third argument is an
// array of function parameters; the null pointer indicates that the
// current set will be used.
   Double_t *params(0);
   return m_func->Integral(a, b, params, eps);
}

/**
 * @file PsfModel.cxx
 * @brief Implementation for fitting a parameterization of the psf.
 * @author J. Chiang
 *
 * $Header$
 */

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include "TBranch.h"

#include "PsfModel.h"

PsfModel::PsfModel(const std::string &outputFile, int maxTrys) 
   : Fitter(outputFile), m_maxTrys(maxTrys) {

   m_func 
      = new TF1("psfModel", 
                "[0]*x*exp(-0.5*pow((x-[1])/[2], 2.))*exp(-sqrt(x/[3]))");

   if (m_tree) {
// Set up the branches in the output tree.
      m_params.resize(4);
      m_tree->Branch("prefactor" , &m_params[0], "prefactor/D");
      m_tree->Branch("mean" , &m_params[1], "mean/D");
      m_tree->Branch("sigma" , &m_params[2], "sigma/D");
      m_tree->Branch("scale" , &m_params[3], "scale/D");
   }
}

PsfModel::~PsfModel() {
}

void PsfModel::applyFit(TH1 * h) {
// These initial parameters were obtained by running Root
// interactively on a sample merit file.
   m_params[0] = 1e4;
   m_params[1] = 30.;
   m_params[2] = 10.;
   m_params[3] = 0.05;
   m_func->SetParameters(&m_params[0]);
   m_func->SetParLimits(1, 0, 100);
   int fitTrys = 0;
// // Do the fit in quiet mode.
//    while (h->Fit("psfModel", "Q") && fitTrys++ < m_maxTrys);
   while (h->Fit("psfModel") && fitTrys++ < m_maxTrys);
}

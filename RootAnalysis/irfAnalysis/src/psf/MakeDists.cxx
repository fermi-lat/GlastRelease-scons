/**
 * @file MakeDists.cxx
 * @brief Interface for creating distributions using Root's Project()
 * method.
 * @author J. Chiang
 *
 * $Header$
 */

#include <vector>
#include <iostream>
#include <sstream>

#include "TF1.h"

#include "PSF.h"
#include "Fitter.h"
#include "MakeDists.h"

MakeDists::MakeDists(const std::string &summaryFile) : IRF() {
   m_summary_filename = output_file_root() + summaryFile;
}

void MakeDists::project(const std::string &branchName, 
                        double xmin, double xmax, int nbins,
                        Fitter * fitter) {

   open_input_file();

// Output file for the histograms.
   TFile histFile(m_summary_filename.c_str(), "recreate");
   std::cout << "Writing histograms to " << m_summary_filename 
             << std::endl;

// Loop over angle ranges.
   for (int i = 0; i < angle_bins; ++i) {
      std::cout << angle_cut(i) << std::endl;
      TCut angle(angle_cut(i));
      
// Loop over energies.
      for (int j = 0; j < energy_bins; ++j) {
         double logecenter = logestart + logedelta*j;
         double ecenter = pow(10, logecenter);
         std::ostringstream title;
         title << "Scaled distribution: "
               << (int)(ecenter+0.5) << " MeV, "
               << angles[i] << "--"
               << angles[i+1] << " degrees";
         TCut energy(energy_cut(j));
         TH1F * h = new TH1F(hist_name(i, j), title.str().c_str(), 
                             nbins, xmin, xmax);
         h->GetXaxis()->SetTitle(branchName.c_str());
         std::cout << "\t" << title.str() << "... ";
         m_tree->Project(h->GetName(), branchName.c_str(),
                         goodEvent && energy && angle);

         double scale = h->Integral();
         if (scale > 0) { 
            h->Scale(1./scale);
         }
         double mean = h->GetMean();
         double rms = h->GetRMS();
         std::cout << "count " << scale << ", "
                   << "mean " << mean << ", "
                   << "rms " << rms << std::endl;

         if (fitter) {
            fitter->applyFit(h);
            fitter->writeFitPars();
         }

// Move to the summary file.
         h->SetDirectory(&histFile);
      }
   }
   histFile.Write();
}

void MakeDists::draw(const std::string &ps_filename, double ymax) {
// Delegate to PSF::draw(...) method.
   PSF psf;
   psf.m_summary_filename = m_summary_filename;
   psf.draw(ps_filename,ymax);
}

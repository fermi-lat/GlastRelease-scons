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
#include "TProfile.h"
#include "TPostScript.h"

#include "Fitter.h"
#include "MakeDists.h"

MakeDists::MakeDists(const std::string &summaryFile, bool makeProfile)
   : PSF(), m_makeProfile(makeProfile) {
   m_summary_filename = output_file_root() + summaryFile;
}

void MakeDists::project(const std::string &branchName, 
                        double xmin, double xmax, int nbins,
                        Fitter * fitter) {

    

   m_nbins = nbins;

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
         title << (int)(ecenter+0.5) << " MeV, "
               << angles[i] << "-"
               << angles[i+1] << " degrees";
         TCut energy(energy_cut(j));
         TH1 * h;
         if (m_makeProfile) {
            h = new TProfile(hist_name(i, j), title.str().c_str(), 
                             nbins, xmin, xmax);
         } else {
            h = new TH1F(hist_name(i, j), title.str().c_str(), 
                         nbins, xmin, xmax);
         }
         h->GetXaxis()->SetTitle(branchName.c_str());
         std::cout << "\t" << title.str() << "... ";
         m_tree->Project(h->GetName(), branchName.c_str(),
                         goodEvent && energy && angle);

         if (!m_makeProfile) {
            double scale = h->Integral();
            if (scale > 0) { 
               h->Scale(1./scale);
            }
            double mean = h->GetMean();
            double rms = h->GetRMS();
            std::cout << "count " << scale << ", "
                      << "mean " << mean << ", "
                      << "rms " << rms;
         }
         std::cout << std::endl;

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

void MakeDists::draw(const std::string &ps_filename,  bool logy, Fitter* fitter) {

    TFile psf_file(summary_filename().c_str() ); // for the histograms
    if( ! psf_file.IsOpen()) throw "could not open psf root file";

    if( fitter!=0) gStyle->SetOptFit(111);

    TCanvas c;
    divideCanvas(c,4,2, std::string("Plots from ")
                 +summary_filename());

    TPostScript ps((output_file_root()+ ps_filename).c_str(), 112);

    for( int i=0; i<angle_bins; ++i){
       ps.NewPage();
       for(int j=0; j<energy_bins; ++j){
          c.cd(j+1);
          gPad->SetRightMargin(0.02);
          gPad->SetTopMargin(0.03);
          if(logy)gPad->SetLogy();
          
          double logecenter=logestart+logedelta*j, ecenter=pow(10, logecenter);
          
          double xmin=j/double(energy_bins), xmax= (j+1)/double(energy_bins);
          
// Amazingly, this code also works for Profile plots.
          TH1F* h = (TH1F*)psf_file.Get(hist_name(i,j));
          if( h==0) { 
             std::cerr << "could not find hist " 
                       << hist_name(i,j) << std::endl;
             return;
          }
          // now add overflow to last bin
          h->SetBinContent(m_nbins, h->GetBinContent(m_nbins)
                           +h->GetBinContent(m_nbins+1));
          h->SetMaximum(m_ymax);
          h->SetMinimum(m_ymin);
          h->SetStats(false);
          h->SetLineColor(i+1);
          
          std::cout << "Processing " << h->GetTitle() << std::endl;

         if (fitter!=0) {
            fitter->applyFit(h);
            fitter->writeFitPars();
            h->SetStats(true);
         }

// Rewrite the title to describe the energy and angle ranges.
          char title[256];
          sprintf(title, " %6d MeV, %2d-%2d degrees", 
                  (int)(ecenter+0.5), angles[i], angles[i+1]);
          h->SetTitle(title);   
          h->GetXaxis()->CenterTitle(true);
          h->Draw(); 
       }
       c.Update();
    }
    ps.Close(); // print ps file,
}

void MakeDists::addCutInfo(const std::string &rootFile, 
                           const std::string &treeName) {
// Assume the file is in the root output directory.
   std::string path = ::getenv("output_file_root");
   TFile f( (path + "/" + rootFile).c_str(), "update" );

   TTree * tree = (TTree*)f.Get(treeName.c_str());
   Int_t nentries = (Int_t)tree->GetEntries();

   Double_t angle_min, angle_max, energy;
   TBranch *angleMin = tree->Branch("angle_min", &angle_min, "angle_min/D");
   TBranch *angleMax = tree->Branch("angle_max", &angle_max, "angle_max/D");
   TBranch *ee = tree->Branch("energy", &energy, "energy/D");

   int indx = 0;
   for (int i = 0; i < angle_bins; ++i) {
      angle_min = angles[i];
      angle_max = angles[i+1];
      for (int j =0; j < energy_bins; ++j, indx++) {
         double logecenter = logestart + logedelta*j;
         energy = pow(10, logecenter);
         angleMin->Fill();
         angleMax->Fill();
         ee->Fill();
      }
   }
   tree->Write("", TObject::kOverwrite);
//   f.Close();
}

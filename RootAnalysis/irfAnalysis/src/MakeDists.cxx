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

   m_branchName = branchName;
   if (m_scalingFunction != "") {
      applyEnergyScaling();
      std::cout << "Applying energy scaling: "
                << m_branchName << std::endl;
   }

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
         h->GetXaxis()->SetTitle(m_branchName.c_str());
         std::cout << "\t" << title.str() << "... ";
         m_tree->Project(h->GetName(), m_branchName.c_str(),
             goodEvent && energy && angle);

         std::cout << "count " << h->Integral() << ", "
             << "mean " <<  h->GetMean()<< ", "
             << "rms " << h->GetRMS() << std::endl;

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

void MakeDists::draw(const std::string &ps_filename, bool logy, Fitter* fitter) {

    TFile psf_file(summary_filename().c_str() ); // for the histograms
    if( ! psf_file.IsOpen()) throw "could not open psf root file";

    if( fitter!=0) gStyle->SetOptFit(111);

    TCanvas c;

    TPostScript ps((output_file_root()+ ps_filename).c_str(), 112);

    for( int i=0; i<angle_bins; ++i){
       ps.NewPage();
        divideCanvas(c,4,2, std::string("Plots from ")
                 +summary_filename());
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

          if (!m_makeProfile) {
            double scale = h->Integral();
            if (scale > 0) { 
                h->Sumw2(); // needed to preserve errors
               h->Scale(1./scale);
            }
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

void MakeDists::setEnergyScaling(std::string scalingFunction,
                                 const std::vector<double> &params) {
   m_scalingFunction = scalingFunction;
   m_params = params;
   replace_substring(scalingFunction, "McEnergy", "x");
   m_energyScale = new TF1("energyScale", scalingFunction.c_str());
   m_energyScale->SetParameters(&params[0]);
}

void MakeDists::addEnergyScaling(const std::string &rootFile,
                                 const std::string &treeName) {

   if (m_energyScale) {
// Assume the file is in the root output directory.
      std::string path = ::getenv("output_file_root");
      TFile f( (path + "/" + rootFile).c_str(), "update" );

// Create and name the tree.
      TTree tree("PsfScale", "Energy scaling for McDirErr");
      Double_t energy, scaleFactor;
      tree.Branch("McEnergy", &energy, "McEnergy/D");
      tree.Branch("ScaleFactor", &scaleFactor, "ScaleFactor/D");

// Should be ok to hard-wire these values...
      double emin = 10;   // in MeV
      double emax = 3e5;  // 300 GeV
      int npts = 100;
      double estep = log(emax/emin)/(npts-1);
      for (int i = 0; i < npts; i++) {
         energy = emin*exp(estep*i);
         scaleFactor = m_energyScale->Eval(energy);
         tree.Fill();
      }
      tree.Write();
   }
}

void MakeDists::applyEnergyScaling() {
   replaceVariables();
   modifyBranchName();
}

void MakeDists::replaceVariables() {
   for (unsigned int i = 0; i < m_params.size(); i++) {
      std::ostringstream pvar, pvalue;
      pvar << "[" << i << "]";
      pvalue << m_params[i];
      replace_substring(m_scalingFunction, pvar.str(), pvalue.str());
   }
}

void MakeDists::modifyBranchName() {
   std::ostringstream newBranchName;
   newBranchName << m_branchName 
                 << "/(" << m_scalingFunction << ")";
   m_branchName = newBranchName.str();
}

void MakeDists::replace_substring(std::string &expression,
                                  const std::string &oldstr,
                                  const std::string &newstr) {
   size_t len = oldstr.size();
   size_t pos = expression.find(oldstr);
   while (pos != std::string::npos) {
      expression.replace(pos, len, newstr);
      pos = expression.find(oldstr);
   }
//   std::cout << expression << std::endl;
}

   

/**
 * @file EffArea.cxx
 * @brief Investigate the phi-dependence of the LAT effective area.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <vector>
#include <sstream>

//#include <TAxis.h>

#include "PSF.h"
#include "TH2.h"

/**
 * @class EffArea
 *
 * @brief Create 2D histograms as a function of instrument theta and
 * phi for a logrithmic grid of energy intervals.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class EffArea : public PSF {

public:

   EffArea(const std::string &rootfile) 
      : PSF(rootfile), m_nphi(10), m_nmu(10) {
// Delete the current rootfile, if necessary.
      system( ("rm -f " + summary_filename()).c_str() );

// Open the merit file (and fill the friend tree if necessary).
      open_input_file();
   }

   void makeHists(const std::string &label, TCut cut) {
      TFile hist_file(summary_filename().c_str(), "update");
      std::string title("phi vs cos_theta");
      TH2F * h2 = new TH2F( (label + " " + title).c_str(), title.c_str(),
                            m_nmu, -1., 0., m_nphi, -M_PI, M_PI );
      m_tree->Project( h2->GetName(), "atan2(McYDir, McXDir):McZDir",
                       goodEvent && cut );
      title = "phi";
      TH1F * h1 = new TH1F( (label + " " + title).c_str(), title.c_str(),
                            m_nphi, -M_PI, M_PI );
      m_tree->Project( h1->GetName(), "atan2(McYDir, McXDir)",
                       goodEvent && cut );
      hist_file.Write();
   }

   static void createGrid(std::vector<double> &x, double xmin, double xmax, 
                          int nx, bool makeLog = false) {
      double xstep;
      if (makeLog) {
         xstep = log(xmax/xmin)/(nx - 1);
      } else {
         xstep = (xmax - xmin)/(nx - 1);
      }
      x.reserve(nx);
      for (int i = 0; i < nx; i++) {
         if (makeLog) {
            x.push_back(xmin*exp(xstep*i));
         } else {
            x.push_back(xstep*i + xmin);
         }
      }
   }

private:

   int m_nphi, m_nmu;

};
      
int main() {

// Create the logrithmic energy grid.
   double emin(20.);
   double emax (2e5);
   int nenergies(6);
   std::vector<double> energies;
   EffArea::createGrid(energies, emin, emax, nenergies, true);

   TCut front("Tkr1FirstLayer < 12");
   TCut back("Tkr1FirstLayer >= 12");

   EffArea myEffArea("EffArea.root");

// Loop over energy bands and create the 2D histograms
   for (int i = 0; i < nenergies-1; i++) {
      std::ostringstream energy_cut;
      energy_cut << "McEnergy >= " << energies[i]
                 << " && McEnergy < " << energies[i+1];
      
      TCut energyCut(energy_cut.str().c_str());

      std::ostringstream frontlabel;
      frontlabel << energies[i] << " <= McEnergy < " << energies[i+1]
                 << ", front";

      std::cout << "Working on " << frontlabel.str() << std::endl;
      myEffArea.makeHists(frontlabel.str(), energyCut && front);

      std::ostringstream backlabel;
      backlabel << energies[i] << " <= McEnergy < " << energies[i+1]
                << ", back";

      std::cout << "Working on " << backlabel.str() << std::endl;
      myEffArea.makeHists(backlabel.str(), energyCut && back);
   }
}

/*  
Create a set of histograms to allow analysis of the energy  response
*/
#include "PSF.h"

#include <fstream>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Aeff : public PSF {
public:
    Aeff(std::string filename): PSF(filename){}

    // override this to project only the summary plots of energy for each angle range
    void project() 
    {
        open_input_file();
        TFile  hist_file(summary_filename().c_str(), "recreate"); // for the histograms
        std::cout << " writing aeff plots to " << summary_filename() << std::endl;
        int nbins=40; double xmin=0, xmax=2.0, ymax=0.15;

        for(int i=0; i<angle_bins; ++i){
            // loop over angle ranges
           char title[256];  sprintf(title,"Energy distirbution: for angles  %2d-%2d degrees", angles[i], angles[i+1]);
            // histogram to show energy distribution for each angle
            TH1F* h = new TH1F(hist_name(i,8), title, (5.5-1.0)/0.125, 1.0, 5.5);
            printf("Projecting angle range %s\n",angle_cut(i));
            TCut angle(angle_cut(i));
            m_tree->Project(h->GetName(), "McLogEnergy", goodEvent&&angle);
            h->SetDirectory(&hist_file);  // move to the summary file
        }
        hist_file.Write();
    }

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(){
    Aeff  all("aeff.root");
    std::string psfile(all.output_file_root()+"aeff.ps");

   if( ! all.fileExists())
        all.project();
    all.drawAeff(psfile+"(", "", "Effective Area vs. energy");

    Aeff cut("cut_aeff.root");
    cut.set_user_cut(TCut("BkVeto==0"));
   if( ! cut.fileExists())
        cut.project();
    cut.drawAeff(psfile+")", "", "Effective Area vs energy with background cuts");

    return 0;
}

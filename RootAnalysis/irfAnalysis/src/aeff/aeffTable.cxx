/** @file  aeffTable.cxx  
    @brief Create a set of histograms to allow analysis of the effective area
    $Header$
*/

#include "PSF.h"
#include "TH2.h"
#include <TAxis.h>


#include <fstream>

// Class AeffTable
//    Fill logE vs cosT Histograms for use by latResponse::AeffRootTable
//
//    From requirements               path to input  merit file and filename
//                                    path to output ( => output_file_root )
//    Parameter
//    filename        output .root file for histos   ( => summary_filename )
//    ps_filename     output .ps
//_________________________________________________________________________________

class AeffTable : public PSF {
public:
    AeffTable(std::string filename, std::string ps_filename): PSF(filename)
        , m_ps_filename(ps_filename)
        , m_binsize(0.125)
        , m_zdir_bins(20)    
	, m_emin(0.016), m_emax(160.0)  
	, m_ngen(4660000)
        , m_target_area(6.0)
    {}


    // project the values of Mc logE and cosT, selecting goodEvents
    //__________________________________________________________________________

    void project() 
    {
      open_input_file(); // open merit file and fill friend tree, if not yet existing

        TFile  hist_file(summary_filename().c_str(), "recreate");
        std::cout << "AeffTable: write aeffTables to " << summary_filename() << std::endl;

        std::cout << "AeffTable: project logE vs cosT Histos " << std::endl;

        // determine normalization factor for Aeff
        double  dcosT       = 1.0/m_zdir_bins;
        double  dlogE       = m_binsize;
        double  ngenPerBin  = m_ngen * (dcosT*dlogE)/log10(m_emax/m_emin);
        double  norm_factor = m_target_area/ngenPerBin ;

        std::cout << "AeffTable: Apply normalization factor for Aeff Tables assuming " 
	    << "\n\t" << m_ngen << " events generated uniformly over "
            << "\n\t cos theta    from -1 to 0"
	    << "\n\t logE with E from "<< m_emin << " to " << m_emax << " GeV"
            << "\n\t dcosT     " << dcosT
            << "\n\t dlogE     " << dlogE 
            << "\n\t Result: Evts generated per bin " << ngenPerBin
            << "\n\t         Normalization factor   " << norm_factor
            << std::endl;


        TCut front("Tkr1FirstLayer<12" );
        TCut back( "Tkr1FirstLayer>=12");

        // Project logE versus cosTheta for front, back and front & back
        char title[256];  
        sprintf(title, "logE vs cosTheta, front Tkr" );
        TH2F *hf = new TH2F( "lectf", title, 
                             static_cast<int>((5.5-1.0)/m_binsize), 
                             1.0, 5.5,   m_zdir_bins, -1.0, 0.0 );
        m_tree->Project( hf->GetName(),"McZDir:McLogEnergy", goodEvent && front );
        std::cout << "Total number of front events: " 
                  << hf->Integral() << std::endl;
        hf->Scale(norm_factor);        

        sprintf(title, "logE vs cosTheta, back  Tkr" );
        TH2F *hb = new TH2F( "lectb", title, 
            static_cast<int>((5.5-1.0)/m_binsize), 
                             1.0, 5.5,   m_zdir_bins, -1.0, 0.0 );
        m_tree->Project( hb->GetName(),"McZDir:McLogEnergy", goodEvent && back );
        std::cout << "Total number of back events: " 
                  << hb->Integral() << std::endl;
        hb->Scale(norm_factor);

        sprintf(title, "logE vs cosTheta, front & back" );
        TH2F *ha = new TH2F( "lecta", title, 
            static_cast<int>((5.5-1.0)/m_binsize), 
                             1.0, 5.5,   m_zdir_bins, -1.0, 0.0 );
        m_tree->Project( ha->GetName(),"McZDir:McLogEnergy", goodEvent );
        std::cout << "Total number of combined events: " 
                  << ha->Integral() << std::endl;
        ha->Scale(norm_factor);

        hist_file.Write();
    }

  //  Draw aeff plots and save to .ps
  //__________________________________________________________________________

    void draw(std::string ps, std::string page_title, std::string hist_title)
    {
        TFile hist_file(summary_filename().c_str() ); // for the histograms
        TCanvas c;
        c.SetFillColor(10);
        divideCanvas(c,2,2,page_title + "plots from "+summary_filename());

        TLegend* leg[3]; // need to define leg objects for each plot
	for (int i=0; i<3; ++i ){
            leg[i] = new TLegend(0.1,0.8, 0.30,0.89);
            leg[i]->SetFillColor(10);
            leg[i]->SetBorderSize(1);
            leg[i]->SetTextSize(0.05);           
            if ( i==0 ) { leg[i]->SetHeader(" Front & Back " ); }
	    if ( i==1 ) { leg[i]->SetHeader(" Front "); }
	    if ( i==2 ) { leg[i]->SetHeader(" Back " ); }
        }

        // draw the histos
	std::string namei[3] = {"a", "f", "b"};     
        char name[5];
	for (int i=0; i<3; ++i ){
	    sprintf( name, "lect%s", namei[i].c_str() );
	    std::cout << "AeffTable: Draw Histogram " << name << std::endl;
 
            TH2F* h =(TH2F*)hist_file.Get( name ) ; 
            if(h==0){
                std::cerr << "AEffTable: could not find "<< name 
                          << " in summary file " << hist_file.GetName() <<std::endl;
                return;
            }

            c.cd(i+2);
            printf("\t  drawing %s\n", h->GetTitle());
            // h->Sumw2(); // needed to preserve errors
            h->SetMaximum(1.1);
            h->SetStats(false);
            //use of SetOptStat(10);
            h->SetTitle(hist_title.c_str());

            h->GetXaxis()->SetTitle("log(Egen/ 1MeV)");
            h->GetYaxis()->SetTitle("cosTheta");
            h->GetZaxis()->SetTitle("Aeff (m^2)");

            h->GetXaxis()->SetTitleOffset(1.5);
            h->GetYaxis()->SetTitleOffset(1.8);

            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);
            h->GetZaxis()->CenterTitle(true);

            //h->SetLineWidth(2);

            h->Draw("surf1");
            leg[i]->Draw();
	}
        c.Print(ps.c_str(), "ps" );
        delete[] leg;  
        
    }


    //__________________________________________________________________________

    void doit(){

        if( ! fileExists()) project();

        std::string psfile( output_file_root()+m_ps_filename );
        draw( psfile, "Effective Area as function of LogE and cosT", "" );

        // if more to add to printed file, use multiple calls to draw with
        //   +"(" and +")" added to psfile for start and end 
    }


private:
    std::string m_ps_filename;
    double m_binsize;
    int    m_zdir_bins;

    // the following values should come from a bookkeeping class
    double m_emin, m_emax;
    int    m_ngen;
    double m_target_area;
};


//_____________________________________________________________________________


int main(){

    AeffTable lect( "aeffTable.root", "aeffTable.ps" );
    lect.doit();
    
    return 0;
}

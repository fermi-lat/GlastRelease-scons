/** @file IRF.cxx
  $Header$
  */
  
#include "IRF.h"
#include <iostream>
#include "TPaveLabel.h"
#include "TDatime.h"
#include "TROOT.h"

int IRF::angles[] = {0, 37, 53, 66, 78}; //cos theta 1, .8., .6, .4, .2


IRF::IRF(std::string summary_root_filename)
: m_file(0), m_tree(0)
, m_summary_filename(output_file_root()+summary_root_filename)
, m_ymin(0)
, m_ymax(1)
,m_user_cut("")
{

    // energy binning: 3 per decade
    logestart=5./3., logedelta=1./3.;
    energy_bins = 8;

    // angle binning : delta cos(theta) = 0.2.
    angle_bins=4;

    // redefine some defaults
   gStyle->SetPadColor(10); // out with the gray
   gStyle->SetCanvasColor(10);
   gStyle->SetTitleColor(10);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0); 

   // this applies it to new hist
   gROOT->ForceStyle();
}

void IRF::open_input_file()
{
    static std::string tree_name("1");

     // load the file and access the TTree
    m_file = new TFile(input_filename().c_str());
    if( ! m_file->IsOpen()){
        std::cerr << "Could not open root file " << input_filename() << std::endl;
        throw "could not open input file";
    }
    std::cout << "Opened root file " << input_filename() << std::endl;
    m_tree = (TTree*)m_file->Get(tree_name.c_str());
    if( m_tree==0) {
        std::cerr << "Did not find tree \"" << tree_name << "\" in the input file" << std::endl;
        throw "did not find the TTree";
    }

    // define basic cuts
    goodCal="CalTotRLn>2&&CalEnergySum>5.0 && IMgoodCalProb>0.2";
    goodPSF ="IMcoreProb>0.2&&IMpsfErrPred<3.0"; // this is Bill's minimal cut
    TCut zdir_cut("Tkr1ZDir<-0.2"); 
    TCut bk_cut("BkVeto==0.0");
    goodEvent=goodCal && goodPSF && zdir_cut;

    if( m_user_cut != "") { 
        goodEvent = goodEvent&&m_user_cut;
    }

    std::cout << "good event definition: " << goodEvent.GetTitle() << std::endl;



}

IRF::~IRF(){
    delete m_file;
}

// copy code from TCanvas::Divide, but leave room for a label at top
void IRF::divideCanvas(TCanvas & c, int nx, int ny, std::string top_title)   {
        int color=10;
        double xmargin=0, ymargin=0, top_margin=0.08; 
        c.SetFillColor(color);
        if (nx <= 0) nx = 1;
        if (ny <= 0) ny = 1;
        std::string temp(top_title+" " + TDatime().AsString());

        TPaveLabel*  title = new TPaveLabel( 0.1,0.95, 0.9,0.99, temp.c_str());
         title->SetFillColor(10);
         title->SetBorderSize(0);
         title->Draw();

        Int_t ix,iy;
        Double_t x1,y1,x2,y2;
        Double_t dy = 1/Double_t(ny);
        Double_t dx = 1/Double_t(nx);
        TPad *pad;
        char *tname = new char [strlen(c.GetName())+6];
        Int_t n = 0;
        if (color == 0) color = c.GetFillColor();
        for (iy=0;iy<ny;iy++) {
            y2 = 1 - iy*dy - ymargin;
            y1 = y2 - dy + 2*ymargin;
            if (y1 < 0) y1 = 0;
            if (y1 > y2) continue;
            for (ix=0;ix<nx;ix++) {
                x1 = ix*dx + xmargin;
                x2 = x1 +dx -2*xmargin;
                if (x1 > x2) continue;
                n++;
                sprintf(tname,"%s_%d",c.GetName(),n);
                pad = new TPad(tname,tname,x1,y1*(1-top_margin),x2,y2*(1-top_margin),color);
                pad->SetNumber(n);
                pad->Draw();
            }
        }
        delete [] tname;
        c.Modified();
    }

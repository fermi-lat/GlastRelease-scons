#include "IRF.h"
#include <iostream>
#include "TPaveLabel.h"
#include "TDatime.h"

int IRF::angles[] = {0, 37, 53, 66, 78}; //cos theta 1, .8., .6, .4, .2


IRF::IRF() :m_file(0), m_tree(0)
{
    // define basic cuts
    goodCal="CalTotRLn>2&&CalEnergySum>5.0 && IMgoodCalProb>0.5";
    goodPSF ="IMcoreProb>0.3&&IMpsfErrPred<3.0"; // this is Bill's minimal cut
    TCut tempFix("Tkr1KalEne > 0.5* McEnergy"); // temporary fix for bad energy estimate
    goodEvent=goodCal&&goodPSF&&tempFix;;  

    std::cout << "Applying global cut " << goodEvent.GetTitle() << std::endl;
    // define root files
    file_root=::getenv("file_root");
    input_filename=file_root+"root_files/fulltup.root";

    // energy binning: 3 per decade
    logestart=5./3., logedelta=1./3.;
    energy_bins = 8;

    // angle binning : delta cos(theta) = 0.2.
    angle_bins=4;

     // load the file and access the TTree
    std::cout << "reading from " << input_filename << std::endl;
    m_file = new TFile(input_filename.c_str());
    if( ! m_file->IsOpen()) throw "could not open input file";
    m_tree = (TTree*)m_file->Get("1");
    if( m_tree==0) throw "did not find the TTree";
    
}

IRF::~IRF(){
    delete m_file;
}

// copy code from TCanvas::Divide
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

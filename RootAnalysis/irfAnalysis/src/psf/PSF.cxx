#include "IRF.h"
#include "TProfile.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/* 
plots to study PSF

*/
class PSF : public IRF {
public:
    PSF();
    void project();
    void draw(std::string ps_filename);
    void drawError(std::string ps_filename);
    void drawAsymmetry(std::string ps_filename);
    void drawAeff(std::string ps_filename);

    bool fileExists(){
        TFile f(summary_filename.c_str());
        return f.IsOpen();
    }

    static double probSum[2]; // for defining quantiles
    // histogram and display
    int nbins; 
    double xmin,xmax, ymax;

};
double PSF::probSum[2]={0.68, 0.95}; // for defining quantiles

PSF::PSF()
: nbins(50), xmin(0), xmax(5.0), ymax(0.16)
{
    summary_filename=file_root+"root_files/2/psf_summary.root";
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF::project() 
{
    TFile  psf_file(summary_filename.c_str(), "recreate"); // for the histograms

    std::cout << " writing hists to " << summary_filename << std::endl;

    for(int i=0; i<angle_bins; ++i){

        // loop over angle ranges
        printf("%s\n",angle_cut(i));
        TCut angle(angle_cut(i));

        for(int j=0; j<energy_bins; ++j){
            // loop over energies
            double logecenter=logestart+logedelta*j, ecenter=pow(10, logecenter);
            char title[256];  sprintf(title, "Scaled PSF:  %6d MeV,  %2d-%2d degrees", (int)(ecenter+0.5), angles[i], angles[i+1]);
            TCut energy(energy_cut(j));
            // create and fill a histogram with the scaled dir error
            TH1F* h = new TH1F(hist_name(i,j), title, nbins, xmin, xmax);
            h->SetLineColor(i+1);
            h->GetXaxis()->SetTitle("Scaled deviation");
            printf("\t%s ... ",title);
            m_tree->Project(h->GetName(), "BestDirErr/PSFscaleFactor",  goodEvent&&energy&&angle);
            double scale = h->Integral();
            double quant[2];
            h->GetQuantiles(2,quant, probSum);
            printf(" count %5.0f 68 %5.2f  95 %5.2f\n", scale, quant[0], quant[1]);
            h->SetDirectory(&psf_file);  // move to the summary file
        }


        // special set of plots to display total, errors and asymmetry  for angular range
        int nebins = 16; double emin=5./3., emax=13./3.;

        // histogram to show energy distribution for each angle
        TH1F* h = new TH1F(hist_name(i,8), "McLogEnergy", nebins, emin, emax);
        printf("projecting hist %s, %s\n", h->GetName(), h->GetTitle());
        m_tree->Project(h->GetName(), "McLogEnergy", goodEvent&&angle);

        char qtitle[256]; sprintf(qtitle, "Profile of error for angles %d-%d; McLogEnergy", angles[i],angles[i+1]);
        TProfile * q = new TProfile(hist_name(i,9), qtitle, nebins, emin, emax);
        printf("projecting profile %s, %s\n", q->GetName(), q->GetTitle());
        m_tree->Project(q->GetName(), "PSFscaleFactor:McLogEnergy", goodEvent&&angle);
        q->SetDirectory(&psf_file);

        char ptitle[256]; sprintf(ptitle, "Profile of error ellipse asymmetry for angles %d-%d; McLogEnergy", angles[i],angles[i+1]); 
        TProfile * p = new TProfile(hist_name(i,10), ptitle, nebins, emin, emax);
        printf("projecting profile %s, %s\n", p->GetName(), p->GetTitle());
        m_tree->Project(p->GetName(), "(Tkr1PhiErr-Tkr1ThetaErr)/(Tkr1ThetaErr+Tkr1PhiErr):McLogEnergy", goodEvent&&angle);
        p->SetDirectory(&psf_file);
    }
    psf_file.Write();
}

 
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF::draw(std::string ps_filename)
{

    TFile psf_file(summary_filename.c_str() ); // for the histograms
    if( ! psf_file.IsOpen()) throw "could not open psf root file";

    TCanvas c("psf","psf"); //, 4,2);
    divideCanvas(c,4,2, std::string("Scaled PSF  plots from ")+summary_filename);

    for(int j=0; j<energy_bins; ++j){
        c.cd(j+1);

        TLegend *leg = new TLegend(0.45,0.7, 0.89,0.89);
        leg->SetTextSize(0.04);
        leg->SetHeader("Angle range   68%   95%");
        double logecenter=logestart+logedelta*j, ecenter=pow(10, logecenter);

        double xmin=j/double(energy_bins), xmax= (j+1)/double(energy_bins);

        for( int i=0; i<angle_bins; ++i){
            TH1F* h = (TH1F*)psf_file.Get(hist_name(i,j));
            if( h==0) { std::cerr << "could not find hist " << hist_name(i,j) << std::endl;
                return;
            }
            // now add overflow to last bin
            h->SetBinContent(nbins, h->GetBinContent(nbins)+h->GetBinContent(nbins+1));
            double quant[2];

            h->GetQuantiles(2,quant, probSum);
            double scale= h->Integral();
            if(scale>0) h->Scale(1/scale);
            h->SetMaximum(ymax);
            h->SetStats(false);
            h->SetLineColor(i+1);

            printf("Drawing ...%s\n", h->GetTitle());
            if(i==0){
                char title[256]; // rewrite the title for the multiple plot
                sprintf(title, " %6d MeV", (int)(ecenter+0.5));
                h->SetTitle( title);   
                h->GetXaxis()->CenterTitle(true);
                h->Draw(); 
            } else h->Draw("same");
            char entry[64]; sprintf(entry," %2d - %2d   %5.1f  %5.1f", angles[i], angles[i+1], quant[0],quant[1]);
            leg->AddEntry( h, entry, "l");
        }
        leg->SetFillColor(10);
        leg->SetBorderSize(1);
        leg->Draw();
    }
    c.Print ( ps_filename.c_str() ); // print ps file,
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF::drawError(std::string ps)
{
    TFile psf_file(summary_filename.c_str() ); // for the histograms
    if( ! psf_file.IsOpen()) { return;}

    TCanvas c;
    c.SetFillColor(10);

    gPad->SetLogy();
    TLegend* leg=new TLegend(0.63,0.7, 0.85,0.89);
    leg->SetHeader("Angle ranges  ");
    for( int i=0; i<angle_bins; ++i){
        TProfile * h = (TProfile*)psf_file.Get(hist_name(i,9));
        if( h==0) { std::cerr << "Could not find hist " << hist_name(i,9) << std::endl;
        return;
        }
        h->SetLineColor(i+1);        
        h->SetLineWidth(2);
        h->SetMinimum(1e-3);
        h->SetStats(false);
        h->SetTitle("Scaled Error vs energy");
        h->GetXaxis()->SetTitle("log(Egen/ 1MeV)");
        h->GetXaxis()->CenterTitle(true);
        h->Draw(i==0? "" : "same");
        char entry[16]; sprintf(entry," %d - %d ", angles[i], angles[i+1] );
        leg->AddEntry( h, entry, "l");
    }
    leg->SetFillColor(10);
    leg->SetBorderSize(1);
    leg->Draw();

    c.Print(ps.c_str());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF::drawAsymmetry(std::string ps)
{
    TFile psf_file(summary_filename.c_str() ); // for the histograms
    if( ! psf_file.IsOpen()) { return;}

    TCanvas c;
     c.SetFillColor(10);
    TLegend* leg=new TLegend(0.13,0.7, 0.35,0.89);
    leg->SetHeader("Angle ranges  ");

    for( int i=0; i<angle_bins; ++i){

        TProfile * h = (TProfile*)psf_file.Get(hist_name(i,10));
        if( h==0) { std::cerr << "Could not find hist " << hist_name(i,10) << std::endl;
            return;
        }
        h->SetMaximum(0.5);
        h->SetLineColor(i+1);
        h->SetLineWidth(2);
        h->SetStats(false);
        h->SetTitle("Asymmetry vs energy");
        h->GetXaxis()->SetTitle("log(Egen/ 1MeV)");
        h->GetXaxis()->CenterTitle(true);

        h->Draw(i==0? "" : "same");
        char entry[16]; sprintf(entry," %d - %d ", angles[i], angles[i+1] );
        leg->AddEntry( h, entry, "l");
    }
    leg->SetFillColor(10);
    leg->SetBorderSize(1);
    leg->Draw();

    c.Print(ps.c_str());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PSF::drawAeff(std::string ps)
{
    TFile hist_file(summary_filename.c_str() ); // for the histograms
    TCanvas c;
     c.SetFillColor(10);

    TLegend* leg=new TLegend(0.13,0.7, 0.30,0.89);
    leg->SetHeader("Angle ranges");
     leg->SetTextSize(0.04);

    for(int i=0; i<4; ++i){
            TH1F* h =(TH1F*)hist_file.Get(hist_name(i,8)) ; 
            if(h==0){
                std::cerr << "could not find "<< hist_name(i,8) << " in summary file " << hist_file.GetName() <<std::endl;
                return;
            }
        printf("Drawing %s\n", h->GetTitle());
        h->SetLineColor(i+1);
        h->SetStats(false);
        h->SetTitle("Energy distribution for angular ranges");
        h->GetXaxis()->SetTitle("log(Egen/ 1MeV)");
        h->GetXaxis()->CenterTitle(true);
        h->SetLineWidth(2);

        if(i==0)h->Draw(); else h->Draw("same");
        char entry[16]; sprintf(entry," %d-%d", angles[i], angles[i+1]);
        leg->AddEntry( h, entry, "l");
    }
    leg->SetFillColor(10);
    leg->SetBorderSize(1);
    leg->Draw();
    c.Print(ps.c_str());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
{
    PSF p;
    if(!p.fileExists())  p.project();
    std::string ps = p.file_root+ "canvas_files/2/psf.ps";
    p.draw(ps+"(");
    p.drawAeff(ps);
    p.drawError(ps);
    p.drawAsymmetry(ps+")");

    std::cout << "done" << std::endl;
    return 0;
}

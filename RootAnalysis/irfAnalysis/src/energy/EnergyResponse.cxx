/*  
Create a set of histograms to allow analysis of the energy  response
*/
#include "IRF.h"


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class EnergyResponse : public IRF {
public:
    EnergyResponse();
    void define();
    void draw(std::string ps);
    void drawEnergy(std::string ps);
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EnergyResponse::EnergyResponse()
{
    m_summary_filename=output_file_root()+"energy.root";
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EnergyResponse::define() 
{

    TFile  hist_file(summary_filename().c_str(), "recreate"); // for the histograms

    std::cout << " writing hists to " << summary_filename() << std::endl;

    int nbins=80; double xmin=0, xmax=2.0, ymax=0.15;
    for(int i=0; i<angle_bins; ++i){
        // loop over angle ranges

        // histogram to show energy distribution for each angle
        TH1F* h = new TH1F(hist_name(i,8), "McLogEnergy", (5.5-1.0)/0.125, 1.0, 5.5);
        printf("Projecting angle range %s\n",angle_cut(i));
        TCut angle(angle_cut(i));
        m_tree->Project(h->GetName(), "McLogEnergy", goodEvent&&angle);

        for(int j=0; j<energy_bins; ++j){
            // loop over energies
            char title[256];
            sprintf(title, "Energy response for %d MeV, angles %d-%d", (int)(eCenter(j)+0.5), angles[i],angles[i+1]);
            TH1F* h = new TH1F(hist_name(i,j), title, nbins, xmin, xmax);
            h->SetLineColor(i+1);
            h->GetXaxis()->SetTitle("Emeas/Egen");
            h->GetXaxis()->CenterTitle(true);

            printf("\t %s ... ",title);
            m_tree->Project(h->GetName(), "EvtEnergySumOpt/McEnergy", goodEvent&&TCut(energy_cut(j))&&angle);
            double scale = h->Integral(), mean= h->GetMean(), rms=h->GetRMS();
            printf(" count %5.0f mean %5.2f rms %5.2f\n", scale, mean, rms);
            if(scale>0){ 
                //h->Sumw2();
                h->Scale(1./scale);
                h->SetMaximum(ymax);
                h->SetStats(false);
            }
            h->SetDirectory(&hist_file);  // move to the summary file
        }
    }
    hist_file.Write();
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EnergyResponse::draw(std::string ps)
{

    TFile hist_file(summary_filename().c_str() ); // for the histograms
    if( ! hist_file.IsOpen()){
        std::cerr <<  "could not open psf root file " << summary_filename() << std::endl;
        return;;
    }

    TCanvas c("energy", "energy");
    divideCanvas(c,4,2, std::string("Energy response  plots from ")+summary_filename());
    for(int j=0; j<energy_bins; ++j){

        c.cd(j+1);
        TLegend *leg = new TLegend(0.55,0.7, 0.89,0.89);
        leg->SetTextSize(0.03);
        leg->SetHeader("Angle range mean rms");

        for( int i=0; i<angle_bins; ++i){
            TH1F* h =(TH1F*)hist_file.Get(hist_name(i,j)) ; 
            if(h==0){
                std::cerr << "could not find "<< hist_name(i,j) << " in summary file " << hist_file.GetName() <<std::endl;
                return;
            }
            printf("Drawing ...%s\n", h->GetTitle());
            char title[256]; // rewrite the title for the multiple plot
            sprintf(title, "%6d MeV", (int)(eCenter(j)+0.5));
            h->SetTitle( title); 
            h->GetXaxis()->CenterTitle(true);
            h->SetLineColor(i+1);

            if(i==0)h->Draw(); else h->Draw("same");
            char entry[64]; sprintf(entry," %d - %d   %5.2f %5.2f", 
                angles[i], angles[i+1],h->GetMean(), h->GetRMS());
            leg->AddEntry( h, entry, "l");
        }
        leg->SetFillColor(10);
        leg->SetBorderSize(1);
        leg->Draw();
    }
    c.Print(ps.c_str());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void EnergyResponse::drawEnergy(std::string ps)
{
    TFile hist_file(summary_filename().c_str() ); // for the histograms
    TCanvas c;

    divideCanvas(c,1,1,"Effective Area");

    TLegend* leg=new TLegend(0.13,0.7, 0.30,0.89);
    leg->SetHeader("Angle ranges");
    leg->SetFillColor(10);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.04);

    // determine normalization factor to Aeff
    int ngen=4.66e6; 
    double anglebin=0.2, logebin=0.125, target_area=6., emax=160., emin=0.016;
    double norm_factor=target_area/ngen/anglebin/(logebin/log10(emax/emin));;
    std::cout << "Applying normailzation factor assuming " << ngen 
        << " generated uniformly over:"
        << "\n\t cos theta from 0 to 1"
        << "\n\t log E with E from "<< emin << " to " << emax << std::endl;

    for(int i=0; i<4; ++i){

        TH1F* h =(TH1F*)hist_file.Get(hist_name(i,8)) ; 
        if(h==0){
            std::cerr << "could not find "<< hist_name(i,8) << " in summary file " << hist_file.GetName() <<std::endl;
            return;
        }
        printf("Drawing %s\n", h->GetTitle());
        h->Scale(norm_factor);
        h->SetLineColor(i+1);
        h->SetStats(false);
        h->SetTitle("Effective Area");
        h->GetXaxis()->SetTitle("log(Egen/ 1MeV)");
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->SetTitle("Aeff (m^2)");
        h->GetYaxis()->CenterTitle(true);

        h->SetLineWidth(2);

        if(i==0)h->Draw(); else h->Draw("same");
        char entry[16]; sprintf(entry," %2d-%d ", angles[i], angles[i+1] );
        leg->AddEntry( h, entry, "l");
    }
    leg->Draw();
    c.Print(ps.c_str() );

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(){
    EnergyResponse er;

    if( !er.fileExists() ) er.define();
    if( !er.fileExists() ) { std::cerr << "cound not open  root file " << er.summary_filename() << std::endl;
        return -1;
    }

    std::string psfile(er.output_file_root()+"energy.ps");
    er.draw(psfile+"(");
    er.drawEnergy(psfile+")");

    printf("done\n");
    return 0;
}

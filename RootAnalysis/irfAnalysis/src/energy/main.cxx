/** @file main.h
$Header$
*/

#include "MakeDists.h"
#include "Fitter.h"
#include "TF1.h"

class EnergyModel : public Fitter {
public:

    EnergyModel(const std::string &outputFile, int maxTrys=5) 
        : Fitter(outputFile), m_maxTrys(maxTrys) 
    {
        m_func 
            = new TF1("myModel",  "[0]*exp(-0.5*pow((x-[1])/[2], 2))");
        m_func->SetParLimits(0, 0.0, 1.0);
        m_func->SetParLimits(1, 0.5, 1.5);
	m_func->SetParLimits(2, 0.0, 0.5);
        // Set up the branches in the output tree.
        const char* names[]={"norm", "mean", "sigma"};
        int npars=sizeof(names)/sizeof(void*);
        m_params.resize(npars);
        for( int i=0; i<npars; ++i){
            m_tree->Branch(names[i], &m_params[i], (std::string(names[i])+"/D").c_str());
        }
    }

    void applyFit(TH1 * h) 
    {
        // These initial parameters were obtained interactively
        double pinit[]={1,1.0, 0.2};
        std::copy(pinit, pinit+sizeof(pinit)/sizeof(double), m_params.begin());
        m_func->SetParameters(&m_params[0]);
       int fitTrys = 0;
        while (h->Fit("myModel") && fitTrys++ < m_maxTrys);
    }

private:
    int m_maxTrys;
};

int main(){
    EnergyModel* myfit = new EnergyModel("energy_fit_parameters.root");

    MakeDists irf("energy_fit.root");

   if ( !irf.fileExists() )
        irf.project("EvtEnergySumOpt/McEnergy", 0, 2, 100 );
    irf.set_ymin(1e-3);
    irf.set_ymax(0.2);
    irf.draw("energy_fit.ps",  true , myfit);
    delete myfit;
    irf.addCutInfo("energy_fit_parameters.root", "fitParams");

   EnergyModel * myfit = new EnergyModel("energy_fit_parameters_thick.root");

    EnergyModel * myfit = new EnergyModel("edisp_thick_parameters.root");
    MakeDists irf_thick("energy_fit_thick.root");
    irf_thick.set_user_cut(TCut("Tkr1FirstLayer>=12.0"));
    if(!irf_thick.fileExists())
       irf_thick.project("EvtEnergySumOpt/McEnergy", 0, 2, 100);
    irf_thick.set_ymin(1e-3);
    irf_thick.set_ymax(0.2);
    irf_thick.draw("energy_fit_thick.ps", true, myfit);
    delete myfit;
    irf_thick.addCutInfo("edisp_thick_parameters.root", "fitParams");

    myfit = new EnergyModel("edisp_thin_parameters.root");
    MakeDists irf_thin("energy_fit_thin.root");
    irf_thin.set_user_cut(TCut("Tkr1FirstLayer<12.0"));
    if(!irf_thin.fileExists())
       irf_thin.project("EvtEnergySumOpt/McEnergy", 0, 2, 100);
    irf_thin.set_ymin(1e-3);
    irf_thin.set_ymax(0.2);
    irf_thin.draw("energy_fit_thin.ps", true, myfit);
    delete myfit;
    irf_thin.addCutInfo("edisp_thin_parameters.root", "fitParams");

    return 0;
}


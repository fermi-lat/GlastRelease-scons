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
            = new TF1("myModel",   "x*[0]*exp(-0.5*pow((x-[1])/[2], 2))");

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
    EnergyModel myfit("energy_fit_parameters.root");

    MakeDists irf("energy_fit.root");

   // if ( !irf.fileExists() )
        irf.project("EvtEnergySumOpt/McEnergy", 0, 2, 100 );
    irf.set_ymin(1e-3);
    irf.set_ymax(1.0);
    irf.draw("energy_fit.ps",  true , &myfit);

    return 0;
}


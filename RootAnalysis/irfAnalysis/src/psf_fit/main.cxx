/** @file main.h
$Header$
*/

#include "MakeDists.h"
#include "Fitter.h"
#include "TF1.h"

class PsfModel : public Fitter {
public:

    PsfModel(const std::string &outputFile, int maxTrys=5) 
        : Fitter(outputFile), m_maxTrys(maxTrys) 
    {
        m_func 
            = new TF1("myModel", 
            "x*( [0]*exp(-0.5*x*x/[1]*[1])/[1] + [2]*exp(-x/[3]) + [4]*exp(-x/[5]) )");

        // Set up the branches in the output tree.
        const char* names[]={"a1", "b1", "a2", "b2", "a3", "b3"};
        int npars=sizeof(names)/sizeof(void*);
        m_params.resize(npars);
        for( int i=0; i<npars; ++i){
            m_tree->Branch(names[i], &m_params[i], (std::string(names[i])+"/D").c_str());
        }
    }
    ~PsfModel(){delete m_func;}

    void applyFit(TH1 * h) 
    {
        // These initial parameters were obtained interactively
        double pinit[]={0.05,0.2, 0.2,1.0, 4e-3,3};
        std::copy(pinit, pinit+sizeof(pinit)/sizeof(double), m_params.begin());
        m_func->SetParameters(&m_params[0]);
       int fitTrys = 0;
        while (h->Fit("myModel") && fitTrys++ < m_maxTrys);
    }

private:
    int m_maxTrys;
};

int main(){

    PsfModel myfit("psf_fit_parameters.root");

    MakeDists psf("psffit2.root");

    if ( !psf.fileExists() )
        psf.project("BestDirErr/(0.05*pow(McEnergy/100,-0.66))", 0, 10, 200 );
    psf.draw("psffit.ps", 0.05, true , &myfit);

    return 0;
}


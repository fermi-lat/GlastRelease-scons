/** @file main.h
$Header$
*/

#include "MakeDists.h"
#include "Fitter.h"
#include "TF1.h"
#include <cmath> // for mod

class PsfModel : public Fitter {
public:

    PsfModel(const std::string &outputFile, int maxTrys=5) 
        : Fitter(outputFile), m_maxTrys(maxTrys) 
    {
        m_func 
            = new TF1("myModel", 
            "x*( [0]*exp(-0.5*x*x/[1]*[1])/[1] + [2]*exp(-x/[3]) )");  //+ [4]*exp(-x/[5]) )");

        // Set up the branches in the output tree.
        const char* names[]={"a1", "b1", "a2", "b2"}; //, "a3", "b3"};
        int npars=sizeof(names)/sizeof(void*);
        m_params.resize(npars);
        for( int i=0; i<npars; ++i){
            m_tree->Branch(names[i], &m_params[i], (std::string(names[i])+"/D").c_str());
        }
         m_func->SetParLimits(0, 0, 1);
         m_func->SetParLimits(1, 0.5, 5);
         m_func->SetParLimits(2, 0, 5);;
         m_func->SetParLimits(3, 0.1, 1);;

    }
    ~PsfModel(){delete m_func;}

    void applyFit(TH1 * h) 
    {
        // These initial parameters were obtained interactively
        double pinit[]={0.01,1.0, 0.2,1.0}; //, 4e-3,3};
        std::copy(pinit, pinit+sizeof(pinit)/sizeof(double), m_params.begin());
        m_func->SetParameters(&m_params[0]);
       int fitTrys = 0;
        while (h->Fit("myModel") && fitTrys++ < m_maxTrys);
    }

private:
    int m_maxTrys;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(){

    PsfModel myfit("psf_fit_parameters.root");

    //-------thin --------------
    MakeDists psf_thin("psf_fit_thin.root");
    psf_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
    if ( !psf_thin.fileExists() )
        psf_thin.project("BestDirErr/(0.05*pow(McEnergy/100,-0.66))", 0, 10, 200 );
    psf_thin.set_ymax(0.1); psf_thin.set_ymin(1e-4);
    psf_thin.draw("psf_fit_thin.ps", true , &myfit);

    //-------thick --------------
    MakeDists psf_thick("psf_fit_thick.root");
    psf_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
    if ( !psf_thick.fileExists() )
        psf_thick.project("BestDirErr/(0.15*pow(McEnergy/100,-0.66))", 0, 10, 200 );
    psf_thick.set_ymax(0.2); psf_thick.set_ymin(1e-4);
    psf_thick.draw("psf_fit_thick.ps", true , &myfit);

    return 0;
}


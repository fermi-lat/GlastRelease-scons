/** @file main.h
$Header$
*/

#include "MakeDists.h"
#include "Fitter.h"
#include "GenericFitter.h"
#include "TF1.h"
#include <cmath> // for mod

class PsfModel : public Fitter {
public:

    PsfModel(const std::string &outputFile, int maxTrys=5) 
        : Fitter(outputFile), m_maxTrys(maxTrys) 
    {
        m_func 
            = new TF1("myModel", 
                "x * ([0] * exp (- 0.5  * x * x /( [1] * [1]))  + [2] * exp (- (x * [3])^[4])   )"
            );
            //"x*( [0]*exp(-0.5*x*x/([1]*[1]))/[1] + [2]*exp(-x/[3]) )");  //+ [4]*exp(-x/[5]) )");

        // Set up the branches in the output tree.
        const char* names[]={"p0", "p1", "p2", "p3", "p4"};;
        int npars=sizeof(names)/sizeof(void*);
        m_params.resize(npars);
        for( int i=0; i<npars; ++i){
            m_tree->Branch(names[i], &m_params[i], (std::string(names[i])+"/D").c_str());
        }
         m_func->SetParLimits(0, 0.1, 1);  // minimum for the gaussian component
         m_func->SetParLimits(1, 0.2, 5);
         m_func->SetParLimits(2, 0, 10);;
         m_func->SetParLimits(3, 0, 10);;
         m_func->SetParLimits(4, 0.5, 1.5);;

    }
    ~PsfModel(){delete m_func;}

    void applyFit(TH1 * h) 
    {
        // These initial parameters were obtained interactively
        double pinit[]={0.01,1.0, 1.0,1.0, 1.0}; //, 4e-3,3};
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

    //-------thin --------------
    Fitter * myfit = new PsfModel("psf_fit_thin_parameters.root");
    MakeDists psf_thin("psf_fit_thin.root");
    psf_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
    if ( !psf_thin.fileExists() )
        psf_thin.project("BestDirErr/PSFscaleFactor", 0, 10, 100 );
    psf_thin.set_ymax(0.2); psf_thin.set_ymin(1e-4);
    psf_thin.draw("psf_fit_thin.ps", true , myfit);
    delete myfit;
    psf_thin.addCutInfo("psf_fit_thin_parameters.root", "fitParams");

    //-------thick --------------
    myfit = new PsfModel("psf_fit_thick_parameters.root");
    MakeDists psf_thick("psf_fit_thick.root");
    psf_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
    if ( !psf_thick.fileExists() )
        psf_thick.project("BestDirErr/PSFscaleFactor", 0, 10, 100 );
    psf_thick.set_ymax(0.5); psf_thick.set_ymin(1e-4);
    psf_thick.draw("psf_fit_thick.ps", true , myfit);
    delete myfit;
    psf_thick.addCutInfo("psf_fit_thick_parameters.root", "fitParams");

    //-------test GenericFitter and setEnergyScaling--------------
    std::string myFunction 
       = "x*( [0]*exp(-0.5*x*x/([1]*[1]))/[1] + [2]*exp(-pow((x/[3]), [4])) )";
    double pinit[] = {0.01, 1.0, 0.2, 1.0, 1.0};
    std::vector<double> fitParams(pinit, pinit+sizeof(pinit)/sizeof(double));
    myfit = new GenericFitter(myFunction, fitParams, 
                              "psf_test_parameters.root");
    double lower[] = {0., 0.5, 0., 0.1, 0.};
    std::vector<double> lowerLims(lower, lower+sizeof(lower)/sizeof(double));
    double upper[] = {1., 5., 5., 1., 3.};
    std::vector<double> upperLims(upper, upper+sizeof(upper)/sizeof(double));
    myfit->setBounds(lowerLims, upperLims);

    MakeDists psf_test("psf_test.root");
    psf_test.set_user_cut(TCut("Tkr1FirstLayer<12"));
    std::vector<double> params;
    params.push_back(0.05);
    params.push_back(-0.66);
    psf_test.setEnergyScaling("[0]*pow(McEnergy/100, [1])", params);
    if ( !psf_test.fileExists() )
       psf_test.project("BestDirErr", 0, 10, 100 );
    psf_test.set_ymax(0.5); psf_test.set_ymin(1e-4);
    psf_test.draw("psf_test.ps", true , myfit);
    delete myfit;
    psf_test.addCutInfo("psf_test_parameters.root", "fitParams");
    psf_test.addEnergyScaling("psf_test_parameters.root", "energyScaling");

    return 0;
}


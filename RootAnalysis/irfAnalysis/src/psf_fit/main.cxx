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

    //-------test thin--------------
    std::string myFunction 
       = "x*( [0]*exp(-0.5*x*x/([1]*[1])) + [2]*exp(-pow((x*[3]), [4])) )";
    double pinit[] = {0.01, 1.0, 0.2, 1.0, 1.0};
    std::vector<double> fitParams(pinit, pinit+sizeof(pinit)/sizeof(double));
    Fitter * myfit = new GenericFitter(myFunction, fitParams, 
                                       "psf_test_thin_parameters.root");
    double lower[] = {0., 0.5, 0., 0.1, 0.};
    std::vector<double> lowerLims(lower, lower+sizeof(lower)/sizeof(double));
    double upper[] = {1., 5., 5., 1., 3.};
    std::vector<double> upperLims(upper, upper+sizeof(upper)/sizeof(double));
    myfit->setBounds(lowerLims, upperLims);

    MakeDists psf_test_thin("psf_test_thin.root");
    psf_test_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
    params[0] = 0.075;
    params[1] = -0.85;
    params[2] = 5e-4;
    std::string energyScaling("[0]*pow(1e-2/(1./McEnergy + [2]), [1])");
    psf_test_thin.setEnergyScaling(energyScaling, params);
    if ( !psf_test_thin.fileExists() ) 
       psf_test_thin.project("BestDirErr", 0, 10, 100 );
    psf_test_thin.set_ymax(0.5); psf_test_thin.set_ymin(1e-4);
    psf_test_thin.draw("psf_test_thin.ps", true , myfit);
    delete myfit;
    psf_test_thin.addCutInfo("psf_test_thin_parameters.root", "fitParams");
    psf_test_thin.addEnergyScaling("psf_test_thin_parameters.root", 
                                   "energyScaling");

    //-------test thick-------------
    myfit = new GenericFitter(myFunction, fitParams, 
                              "psf_test_thick_parameters.root");
    myfit->setBounds(lowerLims, upperLims);

    MakeDists psf_test_thick("psf_test_thick.root");
    psf_test_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
    params[0] = 0.15;
    params[1] = -0.85;
    params[2] = 1e-4;
    psf_test_thick.setEnergyScaling(energyScaling, params);
    if ( !psf_test_thick.fileExists() ) 
       psf_test_thick.project("BestDirErr", 0, 10, 100 );
    psf_test_thick.set_ymax(0.5); psf_test_thick.set_ymin(1e-4);
    psf_test_thick.draw("psf_test_thick.ps", true , myfit);
    delete myfit;
    psf_test_thick.addCutInfo("psf_test_thick_parameters.root", "fitParams");
    psf_test_thick.addEnergyScaling("psf_test_thick_parameters.root", 
                                    "energyScaling");

    //-------energy dispersion-------
    std::string gaussian("[0]/[2]*exp(-0.5*pow((x - [1])/[2], 2.))");
    double gauss_params[] = {1., 1., 0.2};

    std::vector<double> gaussParams(gauss_params, gauss_params+3);
    Fitter * myfit = new GenericFitter(gaussian, gaussParams, 
                                       "edisp_thin_parameters.root");
    MakeDists edisp_thin("edisp_thin.root");
    edisp_thin.set_user_cut(TCut("Tkr1FirstLayer<12"));
    if (!edisp_thin.fileExists()) 
       edisp_thin.project("EvtEnergySumOpt/McEnergy", -0.5, 2.5, 100);
    edisp_thin.set_ymax(0.5);
    edisp_thin.set_ymin(1e-4);
    edisp_thin.draw("edisp_thin.ps", true, myfit);
    delete myfit;
    edisp_thin.addCutInfo("edisp_thin_parameters.root", "fitParams");

    gaussParams[0] = 1.;
    gaussParams[1] = 1.;
    gaussParams[2] = 0.2;
    myfit = new GenericFitter(gaussian, gaussParams, 
                              "edisp_thick_parameters.root");
    MakeDists edisp_thick("edisp_thick.root");
    edisp_thick.set_user_cut(TCut("Tkr1FirstLayer>11"));
    if (!edisp_thick.fileExists()) 
       edisp_thick.project("EvtEnergySumOpt/McEnergy", -0.5, 2.5, 100);
    edisp_thick.set_ymax(0.5);
    edisp_thick.set_ymin(1e-4);
    edisp_thick.draw("edisp_thick.ps", true, myfit);
    delete myfit;
    edisp_thick.addCutInfo("edisp_thick_parameters.root", "fitParams");

    return 0;
}


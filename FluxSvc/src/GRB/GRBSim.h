/**
 * GRBSim: Simulator engine of a GRB source
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Cohen-Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */
#include "FluxSvc/mainpage.h"
//#include <cmath>
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBConstants.h"

#ifndef GRBSIM_H
#define GRBSIM_H 1

class GRBSim 
{
 public:
  GRBSim(char);
  ~GRBSim();
  void Start();
  
  GRBConstants *myParam;
  
  std::vector<double> Energy(){return m_energy;}
  std::vector<double> DeltaE(){return m_energy;}
  double Tmax(){return m_tmax;}
  double Area(){return m_Area;}
  void ComputeFlux(double time);
  std::vector<double> Spectrum(){return m_spectrum;}/// converted in ->photons/s/eV/m^2
  double IFlux(double enmin=cst::enph);  // in eV/(m^2 s) 
  double IRate(double enmin=cst::enph);  // in ph/(m^2 s) 
  double IEnergy(double enmin=cst::enph);// in eV/(m^2) 

  std::pair<float,float> GRBdir(){return _grbdir;}

  double Flux(int en){return m_spectrum[en];}/// converted in ->photons/s/eV/m^2
  double Energy (int en) {return m_energy[en];}
  double tmax(){return m_tmax;}
  float DrawPhotonFromSpectrum(std::vector<double>, float u=0.0, double emin=cst::enph);
  
 private:
  //data member
  std::vector<GRBShell*> theShells;
  std::vector<GRBShock*> theShocks;
  std::vector<double> m_energy,m_de,m_spectrum;
  std::pair<float,float> _grbdir;
  double m_dt;
  double m_tmax;
  double m_ftot;
  double m_phtot;
  double m_Area;
};

#endif



/**
 * GRBSim: Simulator engine of a GRB source
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Coen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */
#include "FluxSvc/mainpage.h"

#include "GRBConstants.h"
#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBUniverse.h"


#ifndef GRBSIM_H
#define GRBSIM_H 1

class GRBSim {

  static const GRBConstants cst;
  
public:
  GRBSim();
  ~GRBSim();
  void Start();
  void TotalProperties(int flag=0);
  void Plot();
  void ComputeFlux(double time);
  std::vector<double> Spectrum(){return m_spectrum;}/// converted in ->photons/s/eV/m^2
  double IFlux(double enmin=cst.enph);  // in eV/(m^2 s) 
  double IRate(double enmin=cst.enph);  // in ph/(m^2 s) 
  double IEnergy(double enmin=cst.enph);// in eV/(m^2) 

  double Flux(int en){return m_spectrum[en];}/// converted in ->photons/s/eV/m^2
  double Energy (int en) {return m_energy[en];}
  double tmax(){return m_tmax;}
  float DrawPhotonFromSpectrum(std::vector<double>,double enmin=cst.enph);
  std::vector<double> m_spectrum;
 private:
  //data member
  std::vector<GRBShell*> theShells;
  std::vector<GRBShock*> theShocks;
  
  GRBUniverse *myUni;
  double m_dt;
  double m_tmax;
  double m_energy[cst.enstep+1];
  double m_de[cst.enstep];
  double m_gammae[cst.enstep+1];
  double m_time[cst.nstep];
  double m_flux[cst.enstep][cst.nstep];
  double m_ftot;
  double m_phtot;
};

#endif



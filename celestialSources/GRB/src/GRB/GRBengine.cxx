#include "GRBShell.h"
#include "GRBShock.h"
#include "GRBConstants.h"
#include "GRBengine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

GRBengine::GRBengine(GRBConstants *myParam)
{
  int engine_type   = 3;
  std::string shell_type = "jet";
  switch(engine_type)
    {
    case 1: // jet, no shell evolution
      {
	std::cout<<" CASE 1: no shell evolution, no collision, physical parameters "<<std::endl;  
	m_duration=0.0;
	double BurstDuration     = myParam->Duration();
	if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
	int NumberOfShocks       = myParam->Nshock();
	double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	double ShockTime         = 0.;
	double ei = myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
	int ns=0;
	
	while (ns < NumberOfShocks)
	  {
	    double  ShockDuration = 0.1;
	    if(NumberOfShocks==1) ShockDuration = BurstDuration;
	    double ComputedShockDuration = 0.0;
	    double gamma_min = myParam->GammaMin();
	    //	    double gamma_max = myParam->GammaMax();
	    double gamma = gamma_min;
	    double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	    
	    GRBShell ShockedMaterial(gamma,mass,myParam->Thickness(),myParam->JetRadius(),myParam->ShellRadius(),shell_type);	  
	    ShockedMaterial.setEint(1.0e-3*ei);
	    GRBShock iShock(ShockedMaterial);
	    
	    ComputedShockDuration = iShock.duration();
	    iShock.setTobs(ShockTime);
	    theShocks.push_back(iShock);		
	    if (ShockTime+ComputedShockDuration>m_duration)
	      {
		m_duration = (ShockTime+ComputedShockDuration < 1.5*BurstDuration) ? 
		  ShockTime+ComputedShockDuration : 1.5*BurstDuration;
	      }
	    ns++;
	    ShockTime += timeBetweenShocks; 
	    timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	  }
      }
      break;
    case 2: // 2 shells that shock, jet, no shell evolution
      {
	std::cout<<"CASE 2: no shell evolution, collision , physical parameters "<<std::endl;
	m_duration=0.0;
	double BurstDuration     = myParam->Duration();
	if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
	int NumberOfShocks       = myParam->Nshock();
	double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	double ShockTime         = 0.;
	double ei = myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
	int ns=0;
	
	while (ns < NumberOfShocks)
	  {
	    double  ShockDuration = 0.1;
	    if(NumberOfShocks==1) ShockDuration = BurstDuration;
	    double ComputedShockDuration = 0.0;
	    double gamma_min = myParam->GammaMin();
	    double gamma_max = myParam->GammaMax();
	    
	    double gamma_front = gamma_min+(gamma_max-gamma_min)*RandFlat::shoot(1.0);
	    double gamma_back  = gamma_front + (gamma_max-gamma_front)*RandFlat::shoot(1.0);
	    double mass_front  = .5*ei/(gamma_front*cst::c2); // Front Shell mass
	    double mass_back   = .5*ei/(gamma_back *cst::c2); // Back  Shell mass
	    GRBShell fShell(gamma_front,mass_front,
			    myParam->Thickness(),
			    myParam->JetRadius(),
			    myParam->ShellRadius(),
			    shell_type);	  
	    GRBShell bShell(gamma_back ,mass_back,
			    myParam->Thickness(),
			    myParam->JetRadius(),
			    myParam->ShellRadius(),
			    shell_type);	  
	    
	    GRBShell ShockedMaterial = fShell + bShell;
	    GRBShock iShock(ShockedMaterial);
	    
	    ComputedShockDuration = iShock.duration();
	    iShock.setTobs(ShockTime);
	    theShocks.push_back(iShock);		
	    
	    if (ShockTime+ComputedShockDuration>m_duration)
	      {
		m_duration = (ShockTime+ComputedShockDuration < 1.5*BurstDuration) ? 
		  ShockTime+ComputedShockDuration : 1.5*BurstDuration;
		
	      }
	    ns++;
	    ShockTime += timeBetweenShocks; 
	    timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	  }
      }
      break;
    case 3: // jet, no shell evolution, obseravbles fixed
      {
	std::cout<<" CASE 3: no shell evolution, no collision, observable fixed "<<std::endl;  
	double rise_time   = myParam->RiseTime(); //sec
	double decay_time  = myParam->DecayTime(); //sec
	double peak_energy = myParam->PeakEnergy(); //MeV

	int NumberOfShocks       = myParam->Nshock();
	double ei = myParam->Etot()/NumberOfShocks; // internal energy of the shocked material
	//////////////////////////////////////////////////
	const double C1=3.3e-11*1.01688;
	const double C2=cst::pi*1.e+10;
	const double C3=2.90851206811141e-14*pow((cst::p-2.)/(cst::p-1.),2.)/(2./(1+cos(myParam->JetAngle()))); // *7.78975e-16*72.*1.01;//7.21787e-15*9.;
	double gamma             = pow(pow(peak_energy,2.)*decay_time/(C3*C3*C2),1./5.);
	double thickness         = rise_time*gamma/C1;
	double jet_radius        = sqrt(ei*gamma*decay_time/(thickness*C2));
	double shell_radius      = 0.5*jet_radius;
	myParam->setGammaMin(gamma);
	myParam->setGammaMax(gamma);
	myParam->setThickness(thickness);
	myParam->setJetRadius(jet_radius);
	myParam->setShellRadius(shell_radius);
	
	m_duration=0.0;
	double BurstDuration     = myParam->Duration();
	if (BurstDuration<=0) BurstDuration=getDurationFromBATSE();
	double timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	double ShockTime         = 0.;
	
	int ns=0;
	
	while (ns < NumberOfShocks)
	  {
	    double ComputedShockDuration = 0.0;
	    double mass = (1.-1.0e-3)*ei/(gamma*cst::c2); // Shell mass
	    //	    double internalEnergy  = ei/
	    GRBShell ShockedMaterial(gamma,mass,thickness,jet_radius,shell_radius,shell_type);	  
	    ShockedMaterial.setEint(1.0e-3*ei);
	    GRBShock iShock(ShockedMaterial);
	    
	    ComputedShockDuration = iShock.duration();
	    iShock.setTobs(ShockTime);
	    theShocks.push_back(iShock);		
	    if (ShockTime+ComputedShockDuration>m_duration)
	      {
		m_duration = ShockTime+ComputedShockDuration;
	      }
	    ns++;
	    ShockTime += timeBetweenShocks; 
	    timeBetweenShocks = 2.*BurstDuration/NumberOfShocks*RandFlat::shoot(1.0);
	    std::cout<<" **************************************************"<<std::endl;
	    std::cout<<" Estimated Rise Time      = "<<C1*thickness/gamma<<std::endl;
	    std::cout<<" Estimated Decay Time     = "<<C2*pow(jet_radius,2.)*thickness/(gamma*ei)<<std::endl;
	    std::cout<<" Estimated Break Energy   = "<<C3*pow(gamma,3.)*sqrt(ei/thickness)/(jet_radius)<<std::endl;
	  }
      }
      break;
    }
  //Glast Direction
  //  m_direction=std::make_pair(((RandFlat::shoot(1.0))*1.4)
  //  			  -0.4,(RandFlat::shoot(1.0))*2*M_PI);
  // Galactic (l,b)
  m_direction = std::make_pair(((RandFlat::shoot(1.0))*360.0)-180.,
			       (RandFlat::shoot(1.0)*180.0)-90.);
  m_distance  = myParam->Distance();
}

double GRBengine::getDurationFromBATSE(char* burst_type)
{
  double dur;
  //////////////////////////////////////////////////
  // It determines if, in case of random selection of the parameters,
  // the burst is long or short...
  if (burst_type!="Short" && burst_type!="Long")
    {
      if (RandFlat::shoot(1.0)<=0.3)
	{burst_type="Short";}
      else
	{burst_type="Long";}
    }  
  if(burst_type=="Short")
    {
      double temp=RandGauss::shoot(-3.65755e-01,5.14102e-01);//GRBConstants::SelectGaussRandom(0.0,2.5);
      dur=pow(10.,temp);
    }
  else
    {
      double temp=RandGauss::shoot(1.46819e+00,4.91505e-01);//GRBConstants::SelectGaussRandom(1.0,3.0);
      dur=pow(10.,temp);
    }
  return dur;
}



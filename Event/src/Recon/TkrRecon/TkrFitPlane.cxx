//----------------------------------------------------------------------
//    
//    Implementation of the Kalman Filter Functions 
//               TkrFitPlane
//
//      Original due to Jose Hernando-Angel circa 1997-1999
//      Re-written to combine both X and Y projections
//      W.B. Atwood, SCIPP/UCSC, Nov. 2001 
//      
//-----------------------------------------------------------------------

#include "Event/Recon/TkrRecon/TkrFitPlane.h"

using namespace Event;
    
void TkrFitPlane::initializeInfo(unsigned int hit, unsigned int tower, unsigned int plane, 
                                 AXIS proj, AXIS nextProj, double z, double energy, 
                                 double radLen, double activeDist) {
    m_IDHit = hit;
    m_IDTower = tower;
    m_IDPlane = plane;
    m_zplane = z;
    m_eneplane = energy;
    m_radLen = radLen;
    m_activeDist = activeDist;
}
    
void TkrFitPlane::initializeHits(const TkrFitHit& meas, const TkrFitHit& pred, const TkrFitHit& fit,
                                 const TkrFitHit& smooth, const TkrFitMatrix& material) {

    m_hitmeas = meas;
    m_hitpred = pred;
    m_hitfit = fit;
    m_hitsmooth = smooth;
    m_Qmaterial = material;
}



void TkrFitPlane::removeHit()
{
	m_IDHit   = 0;
	m_IDTower = 0;
	TkrFitPar pnull(0.,0., 0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));
        TkrFitHit temp(TkrFitHit::MEAS,pnull,covnull);
	setHit(temp);
	clean();
}

void TkrFitPlane::clean()
{
	TkrFitPar pnull(0.,0.,0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));
	setHit(TkrFitHit(TkrFitHit::PRED,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::FIT,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::SMOOTH,pnull,covnull));
}

void TkrFitPlane::clear()
{
	TkrFitPar pnull(0.,0.,0.,0.);
	TkrFitMatrix covnull(HepMatrix(4,4,0));

	m_eneplane = 0.;
	m_IDHit = 0xffffffff;
	m_IDPlane  = -1;
	m_IDTower = -1;	
	m_zplane = 0.;

	setHit(TkrFitHit(TkrFitHit::MEAS,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::PRED,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::FIT,pnull,covnull));
	setHit(TkrFitHit(TkrFitHit::SMOOTH,pnull,covnull));
}

void TkrFitPlane::setHit(const TkrFitHit& hit)
{
    TkrFitHit::TYPE type;
    switch (type=hit.getType()){
    case TkrFitHit::PRED:
	m_hitpred=hit;
	break;
    case TkrFitHit::MEAS:
	m_hitmeas=hit;
	break;
    case TkrFitHit::FIT:
	m_hitfit=hit;
	break;
    case TkrFitHit::SMOOTH:
	m_hitsmooth=hit;
	break;
    case TkrFitHit::UNKNOWN:
        break;
    }   
}

TkrFitHit TkrFitPlane::getHit(TkrFitHit::TYPE type) const
{  
    switch (type){
    case TkrFitHit::PRED:
	return TkrFitHit(m_hitpred);
    case TkrFitHit::MEAS:
	return TkrFitHit(m_hitmeas);
    case TkrFitHit::FIT:
	return TkrFitHit(m_hitfit);
    case TkrFitHit::SMOOTH:
	return TkrFitHit(m_hitsmooth);
    case TkrFitHit::UNKNOWN:
        break;
    } 
    return TkrFitHit();
}

Point TkrFitPlane::getPoint(TkrFitHit::TYPE type)const
{
    TkrFitHit hit=getHit(type);
    return Point(hit.getPar().getXPosition(),
                 hit.getPar().getYPosition(),getZPlane());
}

double TkrFitPlane::getDeltaChiEne(TkrFitHit::TYPE type)const
{
// NOTE: THIS MEMBER FUNCTION IS IDENTICAL TO getDeltaChiSq ?????
//    TkrFitHit hit=getHit(type);
//    double delparX=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
//    double delparY=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
//    double sigma2X=m_hitmeas.getCov().getcovX0X0();
//    double sigma2Y=m_hitmeas.getCov().getcovY0Y0();
    
//    double variance=(delparX*delparX)/sigma2X + (delparY*delparY)/sigma2Y;
//    return variance;

      double variance = getDeltaChiSq(type);
      return variance;
}

void TkrFitPlane::setDeltaEne(double ene)

{       
    double radlen = getRadLen();
 //    Code for e+ & e-: average radiative loss
    double factor = exp(-1.*radlen);
    setEnergy(ene*factor);   

/*
 //  Code for muon testing ~ Bethe-Block dE/dx 
#define MUMASS 105.7
    double mu_sq     = MUMASS*MUMASS; 
    double pb_sq     = ene*ene; //Note: for fitting ene = p*Beta
    double p_sq      = pb_sq*(1.+sqrt(1.+ 4.*mu_sq/pb_sq))/2.;
    double ke        = sqrt(mu_sq+p_sq) - MUMASS; 
    double beta_sq   = pb_sq/p_sq; 
    double d_ke      = radlen*18.3/beta_sq;// const. from wallet card est. 
    double ke_next   = ke - d_ke;
    double e_next    = ke_next + MUMASS;
    double p_next_sq = e_next*e_next - mu_sq; 
    double pB_next   = p_next_sq/e_next;
    setEnergy(pB_next);
  */  
}

double TkrFitPlane::getSigma(TkrFitHit::TYPE type) const
{
// NOTE: THIS MEMBER FUNCTION IS IDENTICAL TO getDeltaChiSq ?????
//    double sigma = 1e6;
//    TkrFitHit hit=getHit(type);
//    double delX=hit.getPar().getXPosition()-hitmeas.getPar().getXPosition();
//    double delY=hit.getPar().getYPosition()-hitmeas.getPar().getYPosition();
//    double sigma2X=hit.getCov().getcovX0X0();
//    double sigma2Y=hit.getCov().getcovY0Y0();
    
//    sigma=(delX*delX)/sigma2X + (delY*delY)/sigma2Y;
//    return sigma;
      double sigma = getDeltaChiSq(type);
      return sigma;
}

double TkrFitPlane::getDeltaChiSq(TkrFitHit::TYPE type) const
{  
    TkrFitHit hit=getHit(type);

    double delpar = 0.;
    double sigma2 = 0.; 
    if(m_projection == TkrCluster::X ) {
        delpar=m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
        sigma2=m_hitmeas.getCov().getcovX0X0() - hit.getCov().getcovX0X0();
    }
    else {
        delpar=m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();
        sigma2=m_hitmeas.getCov().getcovY0Y0() - hit.getCov().getcovY0Y0();
    }

    double chi2 = 1e6;
    if(sigma2 > 0) chi2=(delpar*delpar)/sigma2;
/*
    // Try the full 3D version (see Data Analysis Tech. in HEP by Fruthwirth et al) 
    double rx = m_hitmeas.getPar().getXPosition()-hit.getPar().getXPosition();
    double ry = m_hitmeas.getPar().getYPosition()-hit.getPar().getYPosition();

    // Note we have to to the matrix stuff by hand here since we only meas.
    // co-ordinates and not slopes.. 

    double R11 = m_hitmeas.getCov().getcovX0X0() - hit.getCov().getcovX0X0();
    double R33 = m_hitmeas.getCov().getcovY0Y0() - hit.getCov().getcovY0Y0();
    double R13 = -hit.getCov()(1,3); 

    // Take inverse of R matrix
    double detr = R11*R33 - R13*R13;
    double r11  =  R33/detr;
    double r13  = -R13/detr; 
    double r33  =  R11/detr;

    // Form Chisq. Note: when R13 = 0 (then r13 = 0) this resduce to previous 
    double chi2 = rx*rx*r11 + ry*ry*r33 + 2.*r13*rx*ry;
*/
    return chi2;
}

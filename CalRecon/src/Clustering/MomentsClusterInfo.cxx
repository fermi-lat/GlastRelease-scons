
#include "MomentsClusterInfo.h"

/// Constructor
MomentsClusterInfo::MomentsClusterInfo(const ICalReconSvc* calReconSvc) 
                                     : m_calReconSvc(calReconSvc), m_p0(0.,0.,0.), m_minFrac(0.01)
{
    m_calnLayers = m_calReconSvc->getCalNLayers();

    return;
}

/// This makes CalClusters out of associated CalXtalRecData pointers
Event::CalCluster* MomentsClusterInfo::fillClusterInfo(const XtalDataVec* xTalVec)
{
    // Create an output cluster
    Event::CalCluster* cluster = new Event::CalCluster(0);

    cluster->clear();

    double energy = fillLayerData(xTalVec, cluster);

    fillMomentsData(xTalVec, cluster, energy);

    return cluster;
}

double MomentsClusterInfo::fillLayerData(const XtalDataVec* xTalVec, Event::CalCluster* cluster)
{
    //Initialize local variables
    double ene = 0;                                   // Total energy in this cluster
    Vector pCluster(0.,0.,0.);                        // Cluster position
    std::vector<double> eneLayer(m_calnLayers,0.);    // Energy by layer
    std::vector<Vector> pLayer(m_calnLayers);         // Position by layer
    std::vector<Vector> rmsLayer(m_calnLayers);       // rms by layer

    // Compute barycenter and various moments
    // loop over all crystals in the current cluster
    XtalDataVec::const_iterator xTalIter;
    for(xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                  // crystal energy
        Vector pXtal   = recData->getPosition() - m_p0;         // Vector of crystal position
        int    layer   = (recData->getPackedId()).getLayer();   // layer number

        // update energy of corresponding layer
        eneLayer[layer] += eneXtal;
        
        // update average position of corresponding layer
        Vector ptmp = eneXtal*pXtal;
        pLayer[layer] += ptmp;
        
        // Vector containing squared coordinates, weighted by crystal energy 
        Vector ptmp_sqr(ptmp.x()*pXtal.x(), ptmp.y()*pXtal.y(), ptmp.z()*pXtal.z());
 
        // update quadratic spread, which is proportional to eneXtal;
        // this means, that position error in one crystal
        // is assumed to be 1/sqrt(eneXtal) 
        rmsLayer[layer] += ptmp_sqr;
        
        // update energy sum
        ene  += eneXtal;

        // update cluster position
        pCluster += ptmp;
    }

    // Now take the means
    // if energy sum is not zero - normalize cluster position
    if(ene > 0.) pCluster /= ene; 
 	// if energy is zero - set cluster position to non-physical value
    else pCluster = Vector(-1000., -1000., -1000.);
    
    // loop over calorimeter layers
    int i ;
    for( i = 0; i < m_calnLayers; i++)
    {
        // if energy in the layer is not zero - finalize calculations
        if(eneLayer[i]>0)
        {
            // normalize position in the layer
            pLayer[i] *= (1./eneLayer[i]); 
            
            // normalize quadratic spread in the laye
            rmsLayer[i] *= (1./eneLayer[i]);
            
            // Vector containing the squared average position in each component
            Vector sqrLayer(pLayer[i].x()*pLayer[i].x(),
                            pLayer[i].y()*pLayer[i].y(),
                            pLayer[i].z()*pLayer[i].z());
            
            // the precision of transverse coordinate measurement
            // if there is no fluctuations: 1/sqrt(12) of crystal width
            Vector d;
            double csIWidth = m_calReconSvc->getCalCsIWidth() ;
            if(i%2 == 1) d = Vector(csIWidth*csIWidth/12.,0.,0.);
            else d = Vector(0.,csIWidth*csIWidth/12.,0.);
            
            // subtracting the  squared average position and adding
            // the square of crystal width, divided by 12
            rmsLayer[i] += d-sqrLayer;
        }
            
        
        // if energy in the layer is zero - reset position and spread Vectors
        else 
        {
            pLayer[i]   = m_p0;
            rmsLayer[i] = m_p0;
        }



        // Fill the cluster layer data
        Point layerPos(pLayer[i].x(), pLayer[i].y(), pLayer[i].z());

        // Set the data for the vector
        Event::CalClusterLayerData layerData(eneLayer[i], layerPos, rmsLayer[i]);

        cluster->push_back(layerData);
    }
    // Set energy centroid
	Event::CalParams params(ene, 10*ene, pCluster.x(), pCluster.y(), pCluster.z(), 1.,0.,0.,1.,0.,1.,
                                               0., 0., 1.,   1.,0.,0.,1.,0.,1.);
    cluster->initialize(params, 0., 0.);

    return ene;
}

void MomentsClusterInfo::fillMomentsData(const XtalDataVec* xTalVec, Event::CalCluster* cluster, double energy)
{
    //Initialize local variables
    Vector pCluster = cluster->getPosition();

    // Moments Analysis Here
    // This version lifted directly from code supplied to Bill Atwood by Toby Burnett
    // TU 5/24/2005
	Vector moment(0., 0., 0.); 

    double dircos[3][3];
	dircos[0][1] = 0.; 
    dircos[1][1] = 0.; 
    dircos[2][1] = 0.; 

	if(pCluster.z() > -300.) 
    {
	    double m_xcen     = pCluster.x();
        double m_ycen     = pCluster.y();
        double m_zcen     = pCluster.z();
        int    numXtals   = 0;
		double xprm       = 0.,  yprm = 0.,  zprm = 0., Rsq = 0.,
               Ixx        = 0.,  Iyy  = 0.,  Izz  = 0.,
               Ixy        = 0.,  Ixz  = 0.,  Iyz  = 0.;
        double m_Rsq_mean = 0.;

        for(XtalDataVec::const_iterator xTalIter = xTalVec->begin(); xTalIter != xTalVec->end(); xTalIter++)
        {
            // Construct elements of (symmetric) "Inertia" Tensor:
            // See Goldstein, 1965, Chapter 5 (especially, eqs. 5-6,7,22,26).
            // Analysis easy when translated to energy centroid.
            // get pointer to the reconstructed data for given crystal
		    Event::CalXtalRecData* recData = *xTalIter;
        
            // get reconstructed values
            double eneXtal = recData->getEnergy();                // crystal energy
            if(eneXtal > m_minFrac * energy) 
            {
                numXtals++;
                Vector pXtal   = recData->getPosition() - m_p0;         // Vector of crystal position
                int    layer   = (recData->getPackedId()).getLayer();   // layer number
 
                xprm = pXtal.x() - m_xcen;
                yprm = pXtal.y() - m_ycen;
                zprm = pXtal.z() - m_zcen;

                Rsq = xprm*xprm + yprm*yprm + zprm*zprm;
                m_Rsq_mean += Rsq;
                Ixx += (Rsq - xprm*xprm) * eneXtal;
                Iyy += (Rsq - yprm*yprm) * eneXtal;
                Izz += (Rsq - zprm*zprm) * eneXtal;
                Ixy -= xprm*yprm * eneXtal;
                Ixz -= xprm*zprm * eneXtal;
                Iyz -= yprm*zprm * eneXtal;
            }
        }

        if (numXtals > 0) m_Rsq_mean = m_Rsq_mean / numXtals;

        // Fit Calorimeter Track via principal moments approach.
        double resid = 0.;
        double p, q, r, a, b, rad_test, m, psi;

        // Render determinant of Inertia Tensor into cubic form.
        p = - (Ixx + Iyy + Izz);
        q = Iyy*Izz + Iyy*Ixx + Izz*Ixx - (Ixy*Ixy + Iyz*Iyz + Ixz*Ixz);
        r = - Ixx*Iyy*Izz + Ixx*Iyz*Iyz + Iyy*Ixz*Ixz +
              Izz*Ixy*Ixy - 2.*Ixy*Iyz*Ixz;

        // See CRC's Standard Mathematical Tables (19th edition), pp 103-105.
        // The substitution, y = x - p/3 converts  y^3 + p*y^2 + q*y + r = 0
        // to the form  x^3 + a*x + b = 0 .  Then, if b^2/4 + a^3/27 < 0 ,
        // there will be three real roots -- guaranteed since the Inertia Tensor
        // is symmetric.  A second substitution, x = m*cos(psi) , yields the roots.
        a = (3.*q - p*p)/3.;
        b = (2.*p*p*p - 9.*p*q + 27.*r)/27.;

        rad_test = b*b/4. + a*a*a/27.;

        if ((rad_test < 0.) && (Ixy != 0.) && (Ixz != 0.) && (Iyz != 0.))
        {
            // Construct the roots, which are the principal moments.
            m = 2. * sqrt(-a/3.);
            psi = acos( 3.*b/(a*m) ) / 3.;
            moment[0] = m * cos(psi) - p/3.;
            moment[1] = m * cos(psi + 2.*M_PI/3.) - p/3.;
            moment[2] = m * cos(psi + 4.*M_PI/3.) - p/3.;
        }

        // Construct direction cosines; dircos for middle root is parallel to
        // longest principal axis.
        float A, B, C, D;

        for(int iroot=0; iroot < 3; iroot++) 
        {
            A = Iyz * (Ixx - moment[iroot]) - Ixy*Ixz;
            B = Ixz * (Iyy - moment[iroot]) - Ixy*Iyz;
            C = Ixy * (Izz - moment[iroot]) - Ixz*Iyz;

            D = sqrt( 1. / ( 1./(A*A) + 1./(B*B) + 1./(C*C) ) ) / C;
            dircos[0][iroot] = D * C / A;
            dircos[1][iroot] = D * C / B;
            dircos[2][iroot] = D;
        }
	}

    // Fill CalCluster data
    //Event::CalCluster* cl = new Event::CalCluster(ene, pCluster + p0); 
	// WBA: Note pLayer & rmsLayer are hold overs and of questionable use.  Also the defs
	//      of rms_ were pervert from the original (w.r.t. shower axis) to describe
	//      properties w.r.t. Xtal axis - go figure!  
	Vector caldir = Vector(dircos[0][1], dircos[1][1], dircos[2][1]);

	if(caldir[2] < 0.) caldir = -caldir;

	double rms_long  = (moment[0] + moment[2]) / 2.;
	double rms_trans = moment[1];

    Event::CalParams params(energy, 10*energy, pCluster.x(), pCluster.y(), pCluster.z(), 1.,0.,0.,1.,0.,1.,
                                               caldir.x(),   caldir.y(),   caldir.z(),   1.,0.,0.,1.,0.,1.);

    cluster->initialize(params, rms_long, rms_trans);

    return;
}

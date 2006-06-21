#include <stdio.h>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/Digi/CalDigi.h"
#include "CalUtil/CalDefs.h"
#include "CLHEP/Random/RandGauss.h"
#include "EbfCalData.h"

#include <map>


using idents::CalXtalId;

static StatusCode           get (IGlastDetSvc *detSvc, 
                                 MsgStream       &log,
                                 const char     *name,
                                 double          *val);

/**
 *
 *  @fn unsigned int pack (int xtalId, 
 *                         int rngPos,
 *                         int adcPos,
 *                         int rngNeg,
 *                         int adcNeg)
 *  @brief Internal routine to pack the 5 values describing the response
 *         of a CAL log into a single 32-bit word.
 *
 *  @param xtalId  The column number of the struck log
 *  @param rngPos  The ADC range indicator of the positive log end ADC
 *  @param adcPos  The ADC of the positive log end
 *  @param rngNeg  The ADC range indicator of the negative log end ADC
 *  @param adcNeg  The ADC of the negative log end
 *
**/
static inline unsigned int pack (int xtalId, 
                                 int rngPos,
                                 int adcPos,
                                 int rngNeg,
                                 int adcNeg)
{
 return     ((xtalId & 0xf) << 28)                            | 
          ((((rngNeg & 0x3) << 12) | (adcNeg & 0xfff)) << 14) |
           (((rngPos & 0x3) << 12) | (adcPos & 0xfff));
}

 

/**
 *
 *  @fn StatusCode get (IGlastDetSvc *detSvc,
 *                      MsgStream       &log,
 *                      const char     *name,
 *                      double         scale,
 *                      double          *val)
 *  @brief        Fetches the value of the named parameter from the jobOptions
 *                file and returns it as a value in Mev. 
 *
 *  @param detSvc Pointer to service to retrieve constants from the jobOptions
 *  @param log    The message logging facility
 *  @param name   The name of the parameter to get. All parameters are
 *                read as double precision floats.
 *  @param scale  A scale factor applied to the read value.
 *  @param val    Pointer to receive the returned value
 *
 *  @return       Status 
 *
**/
static StatusCode get (IGlastDetSvc *detSvc,
                       MsgStream       &log,
                       const char     *name,
                       double         scale,
                       double          *val)
{
   StatusCode sc = detSvc->getNumericConstByName (std::string (name), val);


   if (sc.isFailure()) 
   {
       log << MSG::ERROR << "Numeric value for " 
                         <<  name 
                         << " not found" << endreq;
   }
   else
   {

      *val *= scale;
       log << MSG::INFO  << "CAL constant "
                         << name
                         << " = "
                         << *val
                         << endreq;
   }
   
   return sc;
}


   

    

/**
 *
 *  @fn StatusCode EbfCalConstants::initialize (IGlastDetSvc *detSvc, 
 *                                              MsgStream       &log)
 *  @brief  Sucks the parameters controlling the Energy to ADC conversion
 *          and the CAL triggering thresholds from the jobOptions file.
 *
 *  @param  detSvc Services to access the jobOptions file
 *  @param  log    Message logging facility
 *  @return Status, successful iff all the parameters are found.
 *
**/
StatusCode EbfCalConstants::initialize (IGlastDetSvc *detSvc, MsgStream &log)
{
    StatusCode sc;

    /*
     | !!! KLUDGE !!!
     | --------------
     | The TriggerAlg code (from which I swiped this snippet) indicates
     | that the energy data to be readin is in GEV. Experimental results
     | indicate otherwise. I strongly suspect they are in MEV already.
     | This needs to be checked, but I can't find where the numbers are
     | specified.
    */
    static double ToMev = 1.0;
    
  
    sc = get (detSvc, log, "cal.maxResponse0",       ToMev, &m_maxEnergy[0]);
    if (!sc) return sc;
    
    sc = get (detSvc, log, "cal.maxResponse1",       ToMev, &m_maxEnergy[1]);
    if (!sc) return sc;
    
    sc = get (detSvc, log, "cal.maxResponse2",       ToMev, &m_maxEnergy[2]);
    if (!sc) return sc;

    sc = get (detSvc, log, "cal.maxResponse3",       ToMev, &m_maxEnergy[3]);
    if (!sc) return sc;

    sc = get (detSvc, log, "cal.pedestal",             1.0, &m_pedestal);
    if (!sc) return sc;

    sc = get (detSvc, log, "cal.maxAdcValue",          1.0, &m_maxAdc);
    if (!sc) return sc;

    sc = get (detSvc, log, "trigger.LOCALthreshold", ToMev, &m_loThreshold);
    if (!sc) return sc;

    sc = get (detSvc, log, "trigger.HICALthreshold", ToMev, &m_hiThreshold);

    /*
     |
     | !!! KLUDGE !!!
     | --------------
     | This is a wacky conversion. The gain is dependent on the
     | pedestal value (the c.m_maxAdc - c.m_pedestal). Stuck with this
     | until I get some consensus from the off-liners.
     |
    */
    m_adcToMev[0] = m_maxEnergy[0] / (m_maxAdc - m_pedestal);
    m_adcToMev[1] = m_maxEnergy[1] / (m_maxAdc - m_pedestal);
    m_adcToMev[2] = m_maxEnergy[2] / (m_maxAdc - m_pedestal);
    m_adcToMev[3] = m_maxEnergy[3] / (m_maxAdc - m_pedestal);
    /*
      printf ("EbfCalData:initialize: m_loThreshold = %g Mev\n"
      "                       m_hiThreshold = %g Mev\n",
      m_loThreshold,
      m_hiThreshold);
    */
    
    
    return sc;
}




/** 
 * 
 *  @fn    void EbfCalData::initialize ()
 *  @brief Initializes a EbfCalData to handle a new event
 *    
**/
void EbfCalData::initialize ()
{

   StatusCode sc;

   unsigned int      itower;
   EbfCalTowerData   *tower;

   m_lo = m_hi = m_msk = m_TotalEnergy = 0;

   /* 
    |  Loop over all the towers 
    */
   for (itower = 0, tower = m_towers; 
        itower < EbfCalData::NumTowers; 
        itower++,   tower++)
   {
       unsigned int      ilayer;
       EbfCalLayerData   *layer;


       /* Zero the trigger summary bits and layer mask for this tower */
       tower->m_lo  = 0;
       tower->m_hi  = 0;
       tower->m_msk = 0;


       /* Loop over all the layers within this tower */
       for (ilayer = 0, layer = tower->m_layers;
            ilayer < EbfCalTowerData::NumLayers;
            ilayer++,   layer++)      
       {
           layer->m_lo    = 0; 
           layer->m_hi    = 0; 
           layer->m_msk   = 0;
           layer->m_nhits = 0;
       }
   }
  
   return;
}






/**
 *
 *  @fn     void EbfCalData::fill (const Event::CalDigiCol &logs,
 *                                 const EbfCalConstants      &c)
 *  @brief  Fills the hits by tower and layer. It is not assumed that
 *          the hits are sorted within a layer end.
 *
 *  @param  logs  The data for the struck logs on this event
 *  @param     c  The constants used to convert Energy <=> ADC and
 *               the triggering thresholds (in Mev)
 *
**/ 
void EbfCalData::fill (const Event::CalDigiCol &logs,
                       const Event::GltDigi &glt,
                       ICalCalibSvc   *calCalibSvc,
                       const EbfCalConstants      &c)
{
   int print = 0;
 
   initialize ();

   // Protect against an empty CAL record, improbable but does occur 
   if (logs.size ())
   {
       double loThreshold = c.m_loThreshold;
       double hiThreshold = c.m_hiThreshold;
       /*
	 if (print)
	 {
	 printf (" T:L:X Rng      Adc     Energy Rng      Adc     Energy\n"
	 " - - - --- -------- ---------- --- -------- ----------\n");
	 }
       */
       
               

       // Loop over all the CAL hits
       for (Event::CalDigiCol::const_iterator it = logs.begin();
            it != logs.end();
            it++)
       {

           /*
            | 
            | Get the readout mode, either BESTRANGE or ALLRANGE 
            | and the packed CAL ID. The ID contains the tower,
            | layer and column information.
            |
            | !!! KLUDGE !!!
            | --------------
            | When I first implemented this the mode was not correctly
            | supported, so I commented it out. This needs to get fixed.
            |
           */

// Get the  CalDigi
           Event::CalDigi                &digi = **it;

// Find out the trigger mode mode==CalXtalId::BESTRANGE (best range)
//                           mode==CalXtalId::ALLRANGE  (all four ranges)
           idents::CalXtalId::CalTrigMode mode = digi.getMode ();
           m_range4 = (mode==idents::CalXtalId::ALLRANGE) ? true:false;

// For testing set this 4 range readout to true
//            m_range4 = true;
           
//            if(m_range4) printf("EbfCalData: Found 4 Range ReadOut\n");
           
// Get the ID for this digi (tower,layer,log)           
           idents::CalXtalId                id = digi.getPackedId();

// If we have four range read out, Loop over ranges
           int lastRange = (m_range4) ? 4 : 1;
           for(int readoutRange=0; readoutRange<lastRange; readoutRange++) { 

// Get the Readout (request which range)
              const Event::CalDigi::CalXtalReadout 
                            *readout = digi.getXtalReadout(readoutRange);

// Get the Readout (force range 0 for tests)
//              const Event::CalDigi::CalXtalReadout 
//                            *readout = digi.getXtalReadout(0);



// Finally unpack the values for this digi                        
              int           rangePos = readout->getRange(idents::CalXtalId::POS);
              unsigned int    adcPos = readout->getAdc  (idents::CalXtalId::POS);
              int           rangeNeg = readout->getRange(idents::CalXtalId::NEG);
              unsigned int    adcNeg = readout->getAdc  (idents::CalXtalId::NEG);

// Decode the tower, layer and log
              unsigned int   towerId = id.getTower ();
              int            layerId = id.getLayer ();
              unsigned int    xtalId = id.getColumn();

              unsigned int  towerMsk = (1 << towerId);
              unsigned int  layerMsk = (1 << layerId);
              unsigned int   xtalMsk = (1 <<  xtalId);

              EbfCalTowerData *tower = &m_towers[towerId];
              EbfCalLayerData *layer = &tower->m_layers[layerId];
              unsigned int      ihit = layer->m_nhits;

              /*
               | Accumulate the bit masks for the
               |   1. Struck towers within this event
               |   2. Struck layers within this tower
               |   3. Struck logs   within this layer
              */
              m_msk        |= towerMsk;
              tower->m_msk |= layerMsk;
              layer->m_msk |=  xtalMsk;


              /* Add the data for this log to the output stream */
              /* Only count the hits once per log */
              if(readoutRange==0) layer->m_nhits   = ihit + 1;
              
              
/* Kluge Up the 4 range readout  while we wait for GLEAM to do it 
/
/       Assume "full" range is 20 bits
/       Assume the ranges are constructed like the following
/               range 0     bits 0 - 11
/               range 1     bits 3 - 14
/               range 2     bits 6 - 17
/               range 3     bits 9 - 20
/
/       This corresponds to a factor of 8 between the ranges.
/
/       Warning...this is not the way it is actually done...
/       I am just thinking about it this way to get something
/       for testing.  We can't fill in bits that we don't know
/       so just shift around the bits we have...this will artificially
/       give too many zeros, but that all we can do now.
/      
*/
/*
              if(readoutRange > 0) {

// First Positive Side              
                  int newRangePos = rangePos + readoutRange;
                  int newAdcPos = 0;
                  if(newRangePos > 3) {
                      newRangePos -= 4;
                      newAdcPos = (adcPos & 0xfff) << (rangePos-newRangePos)*3;
                  } else {
                      newAdcPos = (adcPos & 0xfff) >> readoutRange*3;
                  }
//                  printf("Pos: %i Old Range %i Old adc 0x%3.3x New Range %i New adc 0x%3.3x\n",readoutRange,rangePos,adcPos,newRangePos,newAdcPos);
                  rangePos = newRangePos;
                  adcPos = newAdcPos;

// Now Negative Side              
                  int newRangeNeg = rangeNeg + readoutRange;
                  int newAdcNeg = 0;
                  if(newRangeNeg > 3) {
                      newRangeNeg -= 4;
                      newAdcNeg = (adcNeg & 0xfff) << (rangeNeg-newRangeNeg)*3;
                  } else {
                      newAdcNeg = (adcNeg & 0xfff) >> readoutRange*3;
                  }
//                  printf("Neg: %i Old Range %i Old adc 0x%3.3x New Range %i New adc 0x%3.3x\n",readoutRange,rangeNeg,adcNeg,newRangeNeg,newAdcNeg);
                  rangeNeg = newRangeNeg;
                  adcNeg = newAdcNeg;                                      
              }
*/              
// Store the hits by packing the information
              layer->m_xtals[xtalId][readoutRange] = pack (xtalId,
                                                           rangePos, adcPos,
                                                           rangeNeg, adcNeg);

              /*
               | CAL TRIGGERING CALCULATION
               | ==========================
               | To do the triggering, must convert the ADC reading along
               | with its range back to MEV, the units that the triggering
               | thresholds are stored in.
               |
               | The CAL contributes two triggering bits from each tower
               | to the GEM, the CAL LO and CAL HI bits. Both are the 
               | straight ORs of all the log ends within the tower. That 
               | is, if any log within a tower is over the LO(HI) threshold,
               | that tower contributes one set bit to the CAL LO(HI)
               | GEM mask.
               |
               | The code below assumes that the hiThreshold > loThreshold.
               | skipping the CAL HI determination if the loThreshold
               | is not satisfied.
               |
               |
               | !!! KLUDGE !!!
               | --------------
               | This really does not belong here. It should be part of a
               | GemDigi class. Unfortunately for now there is no such
               | class, so it is patched in here for the time being.
               |
              */
              double  energyPos = c.convertToMev (rangePos, adcPos);
              double  energyNeg = c.convertToMev (rangeNeg, adcNeg);
	      
//	        if (print)
//	        {
//	        printf (" %1.1x %1.1x %1.1x %3d %8d %10.3f %3d %8d %10.3f \n",
//	        towerId,  layerId, xtalId,
//	        rangePos,  adcPos, energyPos,
//	        rangeNeg,  adcNeg, energyNeg);
//	        }
	      




           if(readoutRange==0) {
           

// Get the information from the GltDigi
              bool fhe_n = glt.getCALHItrigger(idents::CalXtalId(towerId,layerId,xtalId,1));           
              bool fhe_p = glt.getCALHItrigger(idents::CalXtalId(towerId,layerId,xtalId,0));           
//              if(fhe_n || fhe_p)printf("GLT/EBF CAL Hi twr %i layer %i col %i Trigger Neg %i Pos %i\n",towerId,layerId,xtalId,fhe_n,fhe_p);

              bool fle_n = glt.getCALLOtrigger(idents::CalXtalId(towerId,layerId,xtalId,1));           
              bool fle_p = glt.getCALLOtrigger(idents::CalXtalId(towerId,layerId,xtalId,0));           
//              if(fle_n || fle_p)printf("GLT/EBF CAL Lo twr %i layer %i col %i Trigger Neg %i Pos %i\n",towerId,layerId,xtalId,fle_n,fle_p);
  
              
// Store total energy in the event
              m_TotalEnergy += (energyPos+energyNeg)/2;
              
              /* Check if either end is above the CAL LO threshold */ 
              if ( fle_n || fle_p )
              {
                  /* Fill the masks for the low threshold */
                  m_lo        |= towerMsk;
                  tower->m_lo |= layerMsk;
                  layer->m_lo |=  xtalMsk;

//                  printf (" Cal Low Found twr %i lay %i col %i\n",towerId,layerId,xtalId);

               }
                  /* Check if either end is above the CAL HI threshold */ 
               if ( fhe_n || fhe_p )
               {
                   /* Fill the masks for the high threshold */
                   m_hi        |= towerMsk;
                   tower->m_hi |= layerMsk;
                   layer->m_hi |=  xtalMsk;

                   //if (print) printf (" Hi: %8.8x", tower->m_hi);
//                   if(!(fle_n || fle_p ))printf (" Cal High Found with no Cal Lo twr %i lay %i col %i\n",towerId,layerId,xtalId);

               }
              
             } //primary readout range only. 
           } //loop over readout ranges.
           //if (print) printf ("\n");
       } // loop over calDigi entries

// If we are in 4-range readout we should fill in pedestals in channels without hits.
// This is no zero suppression. Seems like CalDig should do this.
      if(m_range4) fillWithPedestals(calCalibSvc);

       
   }


   return;
}


StatusCode EbfCalData::fillWithPedestals(ICalCalibSvc   *calCalibSvc) {

    StatusCode sc;
  
// Loop over towers
   for(int twr=0; twr<16; twr++) {
      for(int lyr=0; lyr<8; lyr++) {
         for(int col=0; col<12; col++) {

              EbfCalTowerData *tower = &m_towers[twr];
              EbfCalLayerData *layer = &tower->m_layers[lyr];
              unsigned int      ihit = layer->m_nhits;
// Check to see if this channel was already filled
              if((layer->m_msk >> col) & 0x1) {
//                 printf("twr %i layer %i col %i already filled \n",twr,lyr,col);
              } else {

// Get pedestal values for each range
                for(int rng = 0; rng<4; rng++) {
                  float ped[2];
                  for(int face=0; face<2; face++) {
//                      idents::CalXtalId rngXtalId(twr,lyr,col,face,rng);
                      CalUtil::RngIdx rngXtalId(twr,lyr,col,face,rng);

                      float tmpPed, tmpPedSig, tmp;
                    //-- RETRIEVE PEDS --//
                       const CalibData::Ped *pedData = calCalibSvc->getPed(rngXtalId);
                       tmpPed = pedData->getAvr();
                       tmpPedSig = pedData->getSig();

                      
                      float rnd = RandGauss::shoot();
                      ped[face] = tmpPed+(tmpPedSig*rnd);
                  }                         
//                  printf("twr %i lyr %i col %i Range %i Ped_N %f Ped_P %f\n",twr,lyr,col,rng,ped[1],ped[0]);

// fill in the values for the pedestals
                  m_msk        |= (1<<twr);
                  tower->m_msk |= (1<<lyr);
                  layer->m_msk |= (1<<col);


                  /* Add the data for this log to the output stream */
                  /* Only count the hits once per log */
                  if(rng==0) layer->m_nhits   = ihit + 1;

// Store the hits by packing the information
                  layer->m_xtals[col][rng] = pack (idents::CalXtalId(twr,lyr,col),
                                                           rng, (int)ped[0],
                                                           rng, (int)ped[1]);


                } 
              }
             
         }
      }
   }

   return StatusCode::SUCCESS;
} 

/*
StatusCode EbfCalData::retrieveCalib(ICalCalibSvc   *calCalibSvc,
                                     unsigned int twr, int lyr, unsigned int col) {
  StatusCode sc;

  using namespace CalDefs;

   float uldThold[2][4];
   
// Loop over faces
  for (FaceNum face; face.isValid(); face++) {


// Note: These thresholds are for the specific ranges
//       fle is assumed to be range 0  (LEX8)
//       fhe is assumed to be range 2  (HEX8)
//  Since we are going to look at a specific (best) range, we need to
//  adjust these for the range that we are interested in.  First get the
//  default values.  
     CalXtalId faceXtalId(twr,lyr,col,face);
     CalibData::ValSig fle,fhe,lac;
     sc = calCalibSvc->getTholdCI(faceXtalId,fle,fhe,lac);
     if (sc.isFailure()) return sc;
     float fleTh = fle.getVal();
     float fheTh = fhe.getVal();
//   lacThresh[face] = lac.getVal();


     for (RngNum rng; rng.isValid(); rng++) {
      CalXtalId rngXtalId(twr,lyr,col,face,rng);

      float tmpPed, tmpPedSig, tmp;
      //-- RETRIEVE PEDS --//
      sc = calCalibSvc->getPed(rngXtalId,
                                 tmpPed,
                                 tmpPedSig,
                                 tmp);                           
      if (sc.isFailure()) return sc;
      m_ped[face][rng] = (float)tmpPed;
      
      //-- RETRIEVE ULD --//
      CalibData::ValSig uld;
      sc = calCalibSvc->getULDCI(rngXtalId,uld);
      if (sc.isFailure()) return sc;
      uldThold[face][rng] = uld.getVal();
//      printf("Face %i Range %i uld %f\n",(int)face,(int)rng,uld.getVal());

// Ok now adjust the thresholds for the range that we have
// Start with the fle range.
      if(rng == LEX1) {
      
        double x8ADC = uldThold[face][0];
        double tmpDAC;
      // 1st convert to dac
        sc = calCalibSvc->evalDAC(CalXtalId(twr, lyr, col, face, LEX8),
                                  x8ADC, tmpDAC);
        if (sc.isFailure()) return sc;
      
        // 2nd convert to next range adc
        double newADC;
        sc = calCalibSvc->evalADC(CalXtalId(twr, lyr, col, face, rng),
                                  tmpDAC, newADC);
        if (sc.isFailure()) return sc;
      
      
        float rat = newADC/x8ADC;
        m_fleThresh[face][rng] = rat*fleTh;
//        printf("New fle found %f for rng %i ( newADC %f  x8ADC %f rat = %f)\n",m_fleThresh[face][rng],(int)rng,newADC,x8ADC,rat);
      } else if(rng==LEX8) {
         m_fleThresh[face][rng] = fleTh;  //LEX8 Threshold
//         printf("Standard fle Threshold %f\n",fleTh);
      } else {
         m_fleThresh[face][rng] = 0x0;  //These high ranges should automatically set low threshold
                                        // This assumes these ranges have saturated.
      }
      

// End with the fhe range.
      if(rng == HEX1) {
      
        double x8ADC = uldThold[face][2];
        double tmpDAC;
      // 1st convert to dac
        sc = calCalibSvc->evalDAC(CalXtalId(twr, lyr, col, face, HEX8),
                                  x8ADC, tmpDAC);
        if (sc.isFailure()) return sc;
      
        // 2nd convert to next range adc
        double newADC;
        sc = calCalibSvc->evalADC(CalXtalId(twr, lyr, col, face, rng),
                                  tmpDAC, newADC);
        if (sc.isFailure()) return sc;
      
      
        float rat = newADC/x8ADC;
        m_fheThresh[face][rng] = rat*fheTh;
//        printf("New fhe found %f for rng %i ( newADC %f  x8ADC %f rat = %f)\n",m_fheThresh[face][rng],(int)rng,newADC,x8ADC,rat);
      } else if (rng==HEX8) {
         m_fheThresh[face][HEX8] = fheTh;  //HEX8 Threshold
//         printf("Standard fhe Threshold %f\n",fheTh);
      } else {
         m_fheThresh[face][rng] = 0xfffffff; // these ranges not used for high threshold
      }

    } //over ranges
   }  // over faces
   
      return StatusCode::SUCCESS;
  
}
*/


/**
 *
 *  @fn     void EbfCalData::fillEncode (const Event::CalDigiCol &logs,
 *                                 const EbfCalConstants      &c)
 *  @brief  Fills the hits by tower and layer. It is not assumed that
 *          the hits are sorted within a layer end.
 *
 *  @param  logs  The data for the struck logs on this event
 *  @param     c  The constants used to convert Energy <=> ADC and
 *               the triggering thresholds (in Mev)
 *
**/ 
void EbfCalData::fillEncode (int encodeFlag, const EbfCalConstants &c, int event)
{
   int print = 0;
 
   initialize ();
   if(encodeFlag==0) return;

   double loThreshold = c.m_loThreshold;
   double hiThreshold = c.m_hiThreshold;      
   
// Loop over each possible log in each layer in each tower
    for (int twr=0; twr<16; twr++){
      for(int layer=0; layer<8; layer++) {
        for(int log=0; log<12; log++) {


// For Small Set (encodeFlag = 2) Only do One Log per tower...but cycle through
	      if( (layer != (event%8) | log != ((event/8)%12) ) && encodeFlag==2 ) continue;

// Loop over ranges
//         m_range4=true; // for testing
         int lastRange = (m_range4) ? 4 : 1;
         for(int readoutRange=0; readoutRange<lastRange; readoutRange++) { 
        

          int           rangePos = readoutRange;
          unsigned int    adcPos = (twr & 0xf)<<8 |  (layer & 0x7)<<4 | (log & 0xf);
        // Make it so special case twr=lay=log=0 does not have a 0 adc count.
          if((twr+layer+log)==0) adcPos |= 1<<7;
          int           rangeNeg = readoutRange;
          unsigned int    adcNeg = (twr & 0xf)<<8 | 1<<7 | (layer & 0x7)<<4 | (log & 0xf);
          printf("ENCODE: Range %i adc+ 0x%8.8x acd- 0x%8.8x\n",readoutRange,adcPos,adcNeg); 

          unsigned int   towerId = twr;
          int            layerId = layer;
          unsigned int    xtalId = log;

// Now have all the information to pack into Ebf Object
          unsigned int  towerMsk = (1 << towerId);
          unsigned int  layerMsk = (1 << layerId);
          unsigned int   xtalMsk = (1 <<  xtalId);

          EbfCalTowerData *tower = &m_towers[towerId];
          EbfCalLayerData *layer = &tower->m_layers[layerId];
          unsigned int      ihit = layer->m_nhits;

        /*
         | Accumulate the bit masks for the
         |   1. Struck towers within this event
         |   2. Struck layers within this tower
         |   3. Struck logs   within this layer
        */
          m_msk        |= towerMsk;
          tower->m_msk |= layerMsk;
          layer->m_msk |=  xtalMsk;


        /* Add the data for this log to the output stream */
          if(readoutRange==0) layer->m_nhits   = ihit + 1;
          layer->m_xtals[xtalId][readoutRange] = pack (xtalId,
                                                     rangePos, adcPos,
                                                     rangeNeg, adcNeg);

        /*
         | CAL TRIGGERING CALCULATION
         | ==========================
         | To do the triggering, must convert the ADC reading along
         | with its range back to MEV, the units that the triggering
         | thresholds are stored in.
         |
         | The CAL contributes two triggering bits from each tower
         | to the GEM, the CAL LO and CAL HI bits. Both are the 
         | straight ORs of all the log ends within the tower. That 
         | is, if any log within a tower is over the LO(HI) threshold,
         | that tower contributes one set bit to the CAL LO(HI)
         | GEM mask.
         |
         | The code below assumes that the hiThreshold > loThreshold.
         | skipping the CAL HI determination if the loThreshold
         | is not satisfied.
         |
         |
         | !!! KLUDGE !!!
         | --------------
         | This really does not belong here. It should be part of a
         | GemDigi class. Unfortunately for now there is no such
         | class, so it is patched in here for the time being.
         |
        */
           double  energyPos = c.convertToMev (rangePos, adcPos);
           double  energyNeg = c.convertToMev (rangeNeg, adcNeg);
	/*
	  if (print)
	  {
	  printf (" %1.1x %1.1x %1.1x %3d %8d %10.3f %3d %8d %10.3f",
	  towerId,  layerId, xtalId,
	  rangePos,  adcPos, energyPos,
	  rangeNeg,  adcNeg, energyNeg);
	  }
	*/


        /* Check if either end is above the CAL LO threshold */ 
          if ( (energyPos > loThreshold) || (energyNeg > loThreshold) )
          {
            /* Fill the masks for the low threshold */
            m_lo        |= towerMsk;
            tower->m_lo |= layerMsk;
            layer->m_lo |=  xtalMsk;

            //if (print) printf (" Lo: %8.8x", tower->m_lo);


            /* Check if either end is above the CAL HI threshold */ 
            if ( (energyPos > hiThreshold) || (energyNeg > hiThreshold) )
            {
                /* Fill the masks for the high threshold */
                m_hi        |= towerMsk;
                tower->m_hi |= layerMsk;
                layer->m_hi |=  xtalMsk;

                //if (print) printf (" Hi: %8.8x", tower->m_hi);
            }
          }

          }//loop over ranges
        } //loop over logs  
      } //loop over layers
    } //loop over towers



   return;
}



/**
 *
 *  @fn         int EbfCalData::format (unsigned int *dst)
 *  @brief      Converts the CAL data into Event Builder Format.
 *
 *  @param  dst The destination array.
 *  @return     Pointer to the next available word in the output array.
 *
 *   It is the responsibility of the caller to ensure that enough room
 *   for a maximum record exists in the output array.  
 *
 *   Since it outputs all the towers as a contigious block, this is 
 *   just a debugging routine. The real formatter must interleave the
 *   CAL and TKR data on a tower-by=tower basis.
 *
**/
unsigned int *EbfCalData::format (unsigned int *dst) const
{
   unsigned int itower;

   /* Loop over each tower */   
   for (itower = 0; itower < EbfCalData::NumTowers; itower++)
   {
       dst = format (dst, itower);
   }

   return dst;
}





/**
 *
 *  @fn unsigned int *EbfCalData::format (unsigned int *dst, int itower) const
 *  @brief  Formats the CAL data for one tower into Event Builder Format
 *
 *  @param    dst The destionation array.
 *  @param itower The tower to format.
 *  @return       Pointer to the next available word in the destination array.
 *
**/
unsigned int *EbfCalData::format (unsigned int *dst, int itower) const
{

   /*
    | This array indicates the order that the layers should be output.
    | and is determined by the GCCC (Calorimeter Cable Controller) on
    | TEM (Tower Electronics Module.
   */
   static const unsigned int Order[EbfCalTowerData::NumLayers] =
   { 
     0, 2, 4, 6,
     1, 3, 5, 7 
   };


   unsigned int            calcnt;
   unsigned int              *beg;
   const EbfCalTowerData   *tower;
   

   /* Check if there is any data in this tower */
   if (( (m_msk & (1 << itower) )== 0) )
   {
       /* No data, 0 the count word and return */
       *dst++ = 0;
       return dst;
   }
        
       
       
   /* Event contains calorimeter data */
   calcnt = 0;               // Zero the hits/layer word 
   beg    = dst++;           // Remember where that hits/layer word goes
   tower  = &m_towers[itower];
   
/* 
   // Loop over the layers in this tower 
   for (unsigned int ilayer = 0; 
        ilayer < EbfCalTowerData::NumLayers;
        ilayer++)
   {

       // Process the layers in the order specified by the 'Order' array 
       int                   olayer = Order[ilayer];


       // Check if there is any data in this layer 
       if ((tower->m_msk & (1 << olayer)) == 0) continue;
       

       // Layer has data, assured that there is at least one log hit 
       const EbfCalLayerData *layer = &tower->m_layers[olayer];
       unsigned int         xtalMsk =  layer->m_msk;


       // Locate the data, keep track of the number of hits this layer 
//       const unsigned int *data = layer->m_xtals;
       int                nhits = 0;
       int               xtalId = 0;
           
             
       // Pack the crystal hits 
       do
       {
           if (xtalMsk & 1)
           {
               nhits += 1;
// Introduce an error so we can see it in the comparison code
//               if(ilayer==3) {
//                  *dst++ = data[xtalId] & 0xffffffa5;
//               } else {   
//                  *dst++  = data[xtalId][0];
                  *dst++  = layer->m_xtals[xtalId][0];
//               }   
           }
       }
       while (xtalId++, xtalMsk >>= 1);

       // 
       // Pack the number of CAL hits one per nibble.
       // This is only necessary when there are hits on a layer, 
       // which is why it is contained within the code that
       // that guarantees there are hits on this layer.
       //
       // Note that the correct layer number is the one that
       // counts which layer currently being output, not the
       // number of the layer that is being output.
       //
       calcnt |= (nhits << (ilayer * 4));
       
   }

*/

   // Loop over groups of two layers (x0,x1),(x2,x3),(y0,y1),(y2,y3). This is
   // important for 4-range readout which groups the information in this fashion  
   for (unsigned int igroup = 0; 
        igroup < EbfCalTowerData::NumLayers/2;
        igroup++){
        
     unsigned int maxRange = (m_range4) ? 4:1;  
     for(unsigned int range = 0; range<maxRange; range++) {   

      for(unsigned int subLayer=0; subLayer<2; subLayer++) {        

       // First which layer does this correspond to
       int ilayer = igroup*2 + subLayer;
       

       // map this choice to the appropriate layer //
       int                   olayer = Order[ilayer];
       //if(itower==1) printf("igrp %i isubL %i ilayer %i  olyr %i rng %i\n",igroup,subLayer,ilayer,olayer,range);

       /* Check if there is any data in this layer */
       if ((tower->m_msk & (1 << olayer)) == 0) continue;
       

       /* Layer has data, assured that there is at least one log hit */
       const EbfCalLayerData *layer = &tower->m_layers[olayer];
       unsigned int         xtalMsk =  layer->m_msk;


       /* Locate the data, keep track of the number of hits this layer */
//       const unsigned int *data = layer->m_xtals;
       int                nhits = 0;
       int               xtalId = 0;
           
             
       /* Pack the crystal hits */
       do
       {
           if (xtalMsk & 1)
           {
               if(range==0) nhits += 1;
// Introduce an error so we can see it in the comparison code
//               if(ilayer==3) {
//                  *dst++ = data[xtalId] & 0xffffffa5;
//               } else {   
//                  *dst++  = data[xtalId][0];
                  *dst++  = layer->m_xtals[xtalId][range];
//               }   
           }
       }
       while (xtalId++, xtalMsk >>= 1);

       /* 
        | Pack the number of CAL hits one per nibble.
        | This is only necessary when there are hits on a layer, 
        | which is why it is contained within the code that
        | that guarantees there are hits on this layer.
        |
        | Note that the correct layer number is the one that
        | counts which layer currently being output, not the
        | number of the layer that is being output.
       */
       calcnt |= (nhits << (ilayer * 4));
      } //subLayer 
     } //4 ranges  
   }//grouped layer


       
   /* Stash the number of CAL hits/layer back in this first word of the dst */
   *beg = calcnt;

   return dst;
}


void EbfCalData::parseInput(unsigned int *contrib, unsigned int itower, unsigned int lcbWords, const EbfCalConstants *calCon)  
{


   static const unsigned int Order[EbfCalTowerData::NumLayers] =
   { 
     0, 2, 4, 6,
     1, 3, 5, 7 
   };

   bool debug = true;
   unsigned int *data = contrib;
   if(debug) {
      printf("Tower %i Contribution:\n",itower);
      for(int i=0; i<lcbWords; i++) {
         for(int j=0; j<4; j++) {
             printf(" 0x%8.8x ",*data++);
         }
         printf("\n");
      }
      printf("\n");
   // reset pointer   
      data = contrib;
   }

// Initialize the EbfCal Class
   if(itower==0) initialize ();
   
   double loThreshold = calCon->m_loThreshold;
   double hiThreshold = calCon->m_hiThreshold;
   
// Peel off header words   
   unsigned int LCB_Header0 = *data++;
   unsigned int LCB_Header1 = *data++;

// Parse out header words   
   unsigned int startbit  = (LCB_Header1 >> 31) & 0x1;
              m_calStrobe = (LCB_Header1 >> 30) & 0x1;
   unsigned int sequence0 = (LCB_Header1 >> 28) & 0x3;
   unsigned int trigAck   = (LCB_Header1 >> 27) & 0x1;
              m_range4    = (LCB_Header1 >> 26) & 0x1;
   unsigned int zeroSup   = (LCB_Header1 >> 25) & 0x1; 
   unsigned int marker    = (LCB_Header1 >> 22) & 0x7;
   unsigned int ErrCont   = (LCB_Header1 >> 21) & 0x1;
   unsigned int DiagCont  = (LCB_Header1 >> 20) & 0x1;
   unsigned int sequence1 = (LCB_Header1 >> 1) & 0x7fff;
   unsigned int ParErr    = (LCB_Header1 & 0x1);

    if(debug) printf("LCB_Header1 0x%8.8x\n",LCB_Header1);
    if(debug) printf("startbit %1i  calStrobe %1i seq %1i trigAck %1i 4Rng %1i ZSup %1i \nMkr %1i ErrCnt %1i DiaCnt %1i ParEr %1i Seq %i\n",
            startbit,m_calStrobe,sequence0,trigAck,m_range4,zeroSup,marker,ErrCont,DiagCont,ParErr,sequence1);
                
// Get Hit words
   unsigned int calHits_Wd = *data++;
   unsigned int calHits[8];
   if(debug) printf("Hits in Cal Layers--");
   for(int i=0; i<8; i++) {
     calHits[i] = (calHits_Wd >> i*4) & 0xf;   
     if(debug) printf(" %i:%i ",i,calHits[i]);
   }
   if(debug) printf("\n");
   
// Strips off hits in each layer
   int nhit=0;
//   int col[12*8];
//   int negRng[12*8];
//   int negADC[12*8];
//   int posRng[12*8];
//   int posADC[12*8];

   int col;
   int negRng;
   int negADC;
   int posRng;
   int posADC;   
   
//   for(int ilayer=0; ilayer<8; ilayer++) {
//      for(int hit=0; hit<calHits[ilayer]; hit++){
// Loop over groups of two layers (x0,x1),(x2,x3),(y0,y1),(y2,y3). This is
// important for 4-range readout which groups the information in this fashion  
   for (unsigned int igroup = 0; 
        igroup < EbfCalTowerData::NumLayers/2;
        igroup++){
        
     unsigned int maxRange = (m_range4) ? 4:1;  
     for(unsigned int range = 0; range<maxRange; range++) {   

      for(unsigned int subLayer=0; subLayer<2; subLayer++) {        
    
// Get number of hits on this layer
        int ilayer = igroup*2 + subLayer;

// Output Layers
        int olayer = Order[ilayer];    

        for(int hit=0; hit<calHits[ilayer]; hit++){
        
// Get the data and parse it
         unsigned int hit_Wd = *data++;
         col    = (hit_Wd & 0xf0000000) >> 28;
         negRng = (hit_Wd & 0x0c000000) >> 26;
         negADC = (hit_Wd & 0x03ffc000) >> 14;
         posRng = (hit_Wd & 0x00003000) >> 12;
         posADC = (hit_Wd & 0x00000fff);
         if(debug) printf("Hit Found on Col %i ilayer %i olayer %i tower %i igroup %i Range %i negRng %i negADC 0x%8.8x posRng %i posADC 0x%8.8x\n",
                              col,ilayer,olayer,itower,igroup,range,negRng,negADC,posRng,posADC);


// Begin to fill the EbfCal Object
         EbfCalTowerData *tower = &m_towers[itower];
         EbfCalLayerData *layer = &tower->m_layers[olayer];
         unsigned int      ihit = layer->m_nhits;

         unsigned int  towerMsk = (1 << itower);
         unsigned int  layerMsk = (1 << olayer);
         unsigned int   xtalMsk = (1 <<  col);
         
         /*
          | Accumulate the bit masks for the
          |   1. Struck towers within this event
          |   2. Struck layers within this tower
          |   3. Struck logs   within this layer
         */
         m_msk        |= towerMsk;
         tower->m_msk |= layerMsk;
         layer->m_msk |=  xtalMsk;
         layer->m_nhits = ihit + 1;
         if(debug) printf("Set the hit masks: tower 0x%8.8x  layer 0x%8.8x   col 0x%8.8x \n",m_msk,tower->m_msk,layer->m_msk);
         
// Store the hits by packing the information
         layer->m_xtals[col][range] = pack (col,
                                                      posRng, posADC,
                                                      negRng, negADC);
     
         if(debug) printf("Packed the hits\n");
// Trigger Information
         double  energyPos = calCon->convertToMev (posRng, posADC);
         double  energyNeg = calCon->convertToMev (negRng, negADC);

         if(range==0) {
// Store total energy in the event
              m_TotalEnergy += (energyPos+energyNeg)/2;
              
              /* Check if either end is above the CAL LO threshold */ 
              if ( (energyPos > loThreshold) || (energyNeg > loThreshold) )
              {
                  /* Fill the masks for the low threshold */
                  m_lo        |= towerMsk;
                  tower->m_lo |= layerMsk;
                  layer->m_lo |=  xtalMsk;

                  //if (print) printf (" Lo: %8.8x", tower->m_lo);


                  /* Check if either end is above the CAL HI threshold */ 
                  if ( (energyPos > hiThreshold) || (energyNeg > hiThreshold) )
                  {
                      /* Fill the masks for the high threshold */
                      m_hi        |= towerMsk;
                      tower->m_hi |= layerMsk;
                      layer->m_hi |=  xtalMsk;

                      //if (print) printf (" Hi: %8.8x", tower->m_hi);
                  }
              }
          } //primary readout range only. 
        if(debug) printf("Set trigger bits m_hi 0x%8.8x  m_lo 0x%8.8x \n",m_hi,m_lo);

        nhit++;


        }//loop over hits
      } //loop over sublayers
     }//loop over ranges
   }//loop over groups   
   
   return;
}



/**
 *
 *  @fn     void print ()
 *  @brief  Diagnostic print routine 
 * 
 *   Routine to print the hits on each crystal of each tower. This
 *   is mainly used to debug the filling routine.
 *
**/
void EbfCalData::print()
{
   unsigned int          itower;
   const EbfCalTowerData *tower;
   

   /* Print the overall CAL LO and CAL HI masks */
   /*
     std::printf ("CAL All    Towers:%8.8x Lo:Hi %8.8x:%8.8x\n", 
     m_msk,
     m_lo,
     m_hi);
   */
   /* Loop over each tower */
   for (itower = 0, tower = m_towers; 
        itower < EbfCalData::NumTowers;
        itower++,   tower++)
   {
       unsigned int          ilayer;
       const EbfCalLayerData *layer;
       
       /*
	 std::printf ("CAL Twr:%1.1x  Layers:%8.8x Lo:Hi %8.8x:%8.8x\n", 
	 itower,
	 tower->m_msk,
	 tower->m_lo,
	 tower->m_hi);
       */

       /* Loop over each layer end of the current tower */ 
       for (ilayer = 0, layer = tower->m_layers;
            ilayer < EbfCalTowerData::NumLayers;
            ilayer++,   layer++)
       {
           unsigned int xtalMsk = layer->m_msk;

           /* Any crystals hit in this layer */       
           if (xtalMsk)
           {
               /* Locate the data */
//               const unsigned int *xtals = layer->m_xtals;
               int                xtalId = 0;
               int              nmargin;
               int                 ncol;
	       
		 nmargin = printf ("       L%1d  Xtals :%8.8x Lo:Hi %8.8x:%8.8x",
		 ilayer,
		 xtalMsk,
		 layer->m_lo,
		 layer->m_hi);
		 ncol    = nmargin;
		 

	       /* Print the crystal hits */
		 do
		 {
		 if (xtalMsk & 1)
		 {
		 if (ncol > 70) 
		 {
		 ncol = std::printf ("\n%*c", nmargin, ' ')-1;
		 }
		 
		 ncol += std::printf (" %8.8x", layer->m_xtals[xtalId][0]);
		 }
		 }
		 while (xtalId++, xtalMsk >>= 1);
		 std::printf ("\n");
           }

       }

   }
       
   return;
}
           







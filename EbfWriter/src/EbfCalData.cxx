#include <stdio.h>
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/Algorithm.h"

#include "Event/Digi/CalDigi.h"
#include "EbfCalData.h"


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
          ((((rngPos & 0x3) << 12) | (adcPos & 0xfff)) << 14) |
           (((rngNeg & 0x3) << 12) | (adcNeg & 0xfff));
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
   unsigned int      itower;
   EbfCalTowerData   *tower;

   m_lo = m_hi = m_msk = 0;

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
           Event::CalDigi                &digi = **it;
//         idents::CalXtalId::CalTrigMode mode = digi.getMode ();
           idents::CalXtalId                id = digi.getPackedId();

           const Event::CalDigi::CalXtalReadout 
                         *readout = digi.getXtalReadout(0);
           int           rangePos = readout->getRange(idents::CalXtalId::POS);
           unsigned int    adcPos = readout->getAdc  (idents::CalXtalId::POS);
           int           rangeNeg = readout->getRange(idents::CalXtalId::NEG);
           unsigned int    adcNeg = readout->getAdc  (idents::CalXtalId::NEG);


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
           layer->m_nhits         = ihit + 1;
           layer->m_xtals[xtalId] = pack (xtalId,
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
            | to the GLT, the CAL LO and CAL HI bits. Both are the 
            | straight ORs of all the log ends within the tower. That 
            | is, if any log within a tower is over the LO(HI) threshold,
            | that tower contributes one set bit to the CAL LO(HI)
            | GLT mask.
            |
            | The code below assumes that the hiThreshold > loThreshold.
            | skipping the CAL HI determination if the loThreshold
            | is not satisfied.
            |
            |
            | !!! KLUDGE !!!
            | --------------
            | This really does not belong here. It should be part of a
            | GltDigi class. Unfortunately for now there is no such
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
               if ( (energyPos > hiThreshold) || (energyPos > hiThreshold) )
               {
                   /* Fill the masks for the high threshold */
                   m_hi        |= towerMsk;
                   tower->m_hi |= layerMsk;
                   layer->m_hi |=  xtalMsk;

                   //if (print) printf (" Hi: %8.8x", tower->m_hi);
               }
           }

           //if (print) printf ("\n");
       }
   }


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
   
 
   /* Loop over the layers in this tower */
   for (unsigned int ilayer = 0; 
        ilayer < EbfCalTowerData::NumLayers;
        ilayer++)
   {

       /* Process the layers in the order specified by the 'Order' array */
       int                   olayer = Order[ilayer];


       /* Check if there is any data in this layer */
       if ((tower->m_msk & (1 << olayer)) == 0) continue;
       

       /* Layer has data, assured that there is at least one log hit */
       const EbfCalLayerData *layer = &tower->m_layers[olayer];
       unsigned int         xtalMsk =  layer->m_msk;


       /* Locate the data, keep track of the number of hits this layer */
       const unsigned int *data = layer->m_xtals;
       int                nhits = 0;
       int               xtalId = 0;
           
             
       /* Pack the crystal hits */
       do
       {
           if (xtalMsk & 1)
           {
               nhits += 1;
               *dst++  = data[xtalId];
           }
       }
       while (xtalId++, xtalMsk >>= 1);

       /* 
        | Pack the number of CAL hits one per nibble.
        | This is only necessary when there are hits on a layer, 
        | which is why it is contained within the code that
        | that guarantees there are hits on this layer.
        |
        | Note that the order of the packing is from the most
        | significant nibble for the first layer out to least
        | significant nibble for the last layer out.
        |
        | Note that the correct layer number is the one that
        | counts which layer currently being output, not the
        | number of the layer that is being output.
       */
       calcnt |= (nhits << (28 - ilayer * 4));
       
   }

       
   /* Stash the number of CAL hits/layer back in this first word of the dst */
   *beg = calcnt;

   return dst;
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
               const unsigned int *xtals = layer->m_xtals;
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
		 
		 ncol += std::printf (" %8.8x", xtals[xtalId]);
		 }
		 }
		 while (xtalId++, xtalMsk >>= 1);
		 std::printf ("\n");
           }

       }

   }
       
   return;
}
           







#include "stdio.h"

#include "Event/Digi/CalDigi.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "EbfCalData.h"

using idents::CalXtalId;

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
 *  @fn    void EbfCalData::initialize ()
 *  @brief Initializes a EbfCalData to handle a new event
 *    
**/
void EbfCalData::initialize ()
{

   StatusCode sc;

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

           // fill all crystal data w/ 0's
           for (unsigned short nLog = 0; nLog < EbfCalLayerData::NumLogs; nLog++) 
             for (unsigned short nReadout = 0; nReadout < 4; nReadout++)
               layer->m_xtals[nLog][nReadout] = 0;
       }
   }
  
   return;
}


/**
 *
 *  @fn     void EbfCalData::fill (const Event::CalDigiCol &calDigiCol)
 *  @brief  Fills the hits by tower and layer. It is not assumed that
 *          the hits are sorted within a layer end.
 *
 *  @param  calDigiCol  The data for the struck calDigiCol on this event
 *
**/ 
StatusCode EbfCalData::fill (const Event::CalDigiCol &calDigiCol,
                             ICalTrigTool &calTrigTool)
{
   int print = 0;
   
   initialize ();

   /*
     if (print)
     {
     printf (" T:L:X Rng      Adc     Energy Rng      Adc     Energy\n"
     " - - - --- -------- ---------- --- -------- ----------\n");
     }
   */
   
   
   
   // Loop over all the CAL hits
   for (Event::CalDigiCol::const_iterator it = calDigiCol.begin();
        it != calDigiCol.end();
        it++)
     {

       /*
         | 
         | Get the readout mode, either BESTRANGE or ALLRANGE 
         | and the packed CAL ID. The ID contains the tower,
         | layer and column information.
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
           const idents::CalXtalId                id = digi.getPackedId();

// Decode the tower, layer and log
           const unsigned int   towerId = id.getTower ();
           const int            layerId = id.getLayer ();
           const unsigned int    xtalId = id.getColumn();

           const unsigned int  towerMsk = (1 << towerId);
           const unsigned int  layerMsk = (1 << layerId);
           const unsigned int   xtalMsk = (1 <<  xtalId);
           
           EbfCalTowerData *tower = &m_towers[towerId];
           EbfCalLayerData *layer = &tower->m_layers[layerId];
           unsigned int      ihit = layer->m_nhits;
           
           /* Add the data for this log to the output stream */
           /* Only count the hits once per log */
           layer->m_nhits   = ihit + 1;

           /// all 4 cal trigger bits for single crystal
           using namespace CalUtil;
           CalArray<XtalDiode, bool> calTriggerBits;
           const XtalIdx xtalIdx(towerId, layerId, xtalId);
           for (XtalDiode xDiode;
                xDiode.isValid();
                xDiode++) {
             const DiodeIdx diodeIdx(xtalIdx, xDiode);
             
             if (calTrigTool.getTriggerBit(diodeIdx, calTriggerBits[xDiode]).isFailure())
               return StatusCode::FAILURE;
           }

           //              if(fle_n || fle_p)printf("GLT/EBF CAL Lo twr %i layer %i col %i Trigger Neg %i Pos %i\n",towerId,layerId,xtalId,fle_n,fle_p);
           
              
           /* Check if either end is above the CAL LO threshold */ 
           if ( calTriggerBits[XtalDiode(POS_FACE, LRG_DIODE)] || calTriggerBits[XtalDiode(NEG_FACE, LRG_DIODE)] )
             {
               /* Fill the masks for the low threshold */
               m_lo        |= towerMsk;
               tower->m_lo |= layerMsk;
               layer->m_lo |=  xtalMsk;
               
               //                  printf (" Cal Low Found twr %i lay %i col %i\n",towerId,layerId,xtalId);
               
             }
           /* Check if either end is above the CAL HI threshold */ 
           if ( calTriggerBits[XtalDiode(POS_FACE, SM_DIODE)] || calTriggerBits[XtalDiode(NEG_FACE, SM_DIODE)] )
             {
               /* Fill the masks for the high threshold */
               m_hi        |= towerMsk;
               tower->m_hi |= layerMsk;
               layer->m_hi |=  xtalMsk;
               
               //if (print) printf (" Hi: %8.8x", tower->m_hi);
               //                   if(!(fle_n || fle_p ))printf (" Cal High Found with no Cal Lo twr %i lay %i col %i\n",towerId,layerId,xtalId);
               
             }

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
              const int           rangePos = readout->getRange(idents::CalXtalId::POS);
              const unsigned int    adcPos = readout->getAdc  (idents::CalXtalId::POS);
              const int           rangeNeg = readout->getRange(idents::CalXtalId::NEG);
              const unsigned int    adcNeg = readout->getAdc  (idents::CalXtalId::NEG);


              /*
               | Accumulate the bit masks for the
               |   1. Struck towers within this event
               |   2. Struck layers within this tower
               |   3. Struck calDigiCol   within this layer
              */
              m_msk        |= towerMsk;
              tower->m_msk |= layerMsk;
              layer->m_msk |=  xtalMsk;

// Store the hits by packing the information
              layer->m_xtals[xtalId][readoutRange] = pack (xtalId,
                                                           rangePos, adcPos,
                                                           rangeNeg, adcNeg);



           } //loop over readout ranges.
           //if (print) printf ("\n");
       } // loop over calDigi entries

       return StatusCode::SUCCESS;
}

void EbfCalData::fillEncode (const int encodeFlag, 
                             const int event,
                             const bool leTriggerBit,
                             const bool heTriggerBit
                             )
{
   int print = 0;
 
   initialize ();
   if(encodeFlag==0) return;

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
         |   3. Struck calDigiCol   within this layer
        */
          m_msk        |= towerMsk;
          tower->m_msk |= layerMsk;
          layer->m_msk |=  xtalMsk;


        /* Add the data for this log to the output stream */
          if(readoutRange==0) layer->m_nhits   = ihit + 1;
          layer->m_xtals[xtalId][readoutRange] = pack (xtalId,
                                                     rangePos, adcPos,
                                                     rangeNeg, adcNeg);

          const bool fhe_n = heTriggerBit;
          const bool fhe_p = heTriggerBit;
           //              if(fhe_n || fhe_p)printf("GLT/EBF CAL Hi twr %i layer %i col %i Trigger Neg %i Pos %i\n",towerId,layerId,xtalId,fhe_n,fhe_p);

          const bool fle_n = leTriggerBit;
          const bool fle_p = leTriggerBit;
           //              if(fle_n || fle_p)printf("GLT/EBF CAL Lo twr %i layer %i col %i Trigger Neg %i Pos %i\n",towerId,layerId,xtalId,fle_n,fle_p);
          
          
          /* Check if either end is above the CAL LO threshold */ 
          if (fle_p || fle_n )
            {
              /* Fill the masks for the low threshold */
              m_lo        |= towerMsk;
              tower->m_lo |= layerMsk;
              layer->m_lo |=  xtalMsk;
              
              //if (print) printf (" Lo: %8.8x", tower->m_lo);
            }
          
          
          /* Check if either end is above the CAL HI threshold */ 
          if (fhe_p || fhe_n)
            {
              /* Fill the masks for the high threshold */
              m_hi        |= towerMsk;
              tower->m_hi |= layerMsk;
              layer->m_hi |=  xtalMsk;
              
              //if (print) printf (" Hi: %8.8x", tower->m_hi);
            }
          
          }//loop over ranges
        } //loop over calDigiCol  
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
           







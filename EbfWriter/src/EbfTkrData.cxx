#include <stdio.h>
#include "GaudiKernel/Algorithm.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "EbfTkrData.h"
#include "Event/Digi/TkrDigi.h"




static inline unsigned int 
     completeTkrTrg (const EbfTkrTowerData *tkr,
                     unsigned int       ntowers);

static int 
        packAccepts (unsigned int                                      *dst,
                     const unsigned int accepts[EbfTkrTowerData::NumCables]);


static int packLayerHits (unsigned int                  *dst, 
                          int                           boff, 
                          const unsigned short int *stripIds,
                          int                          nhits,
                          int                      cableHits);

static int packTots      (unsigned int                  *dst,
                          int                           boff,
                          int                          ntots,
                          const int                    *tots);





/**
 *
 *  @fn    EbfTkrData::initialize ()
 *  @brief Initializes a EbfTkrData class to handle a new event.
 *
**/
void EbfTkrData::initialize ()
{
   unsigned int    itower;
   EbfTkrTowerData *tower;
 
   m_trigger = 0;
   

   /*
    | Zero the hit maps and number of hits on each layer end of each tower
   */
   for (itower = 0, tower = m_towers; 
        itower < sizeof (m_towers) / sizeof (*m_towers); 
        itower++, tower++)
   {
       /* Zero the bit hit maps for each tower */
       tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::LO] = 0;
       tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::HI] = 0;
       tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::LO] = 0;
       tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::HI] = 0;


       /* Zero the number of hits on each layer end */
       for (unsigned int iends = 0; 
            iends < EbfTkrTowerData::NumLayerEnds;
            iends++)
       {
          tower->m_nhits[iends] = 0;
       }
       
   }
   return;
}






/**
 *
 *  @fn     EbfTkrData::fill (const SiData *tkrDigiCol &tkr)
 *  @brief  Fills the hits by tower and layer end. It is assumed 
 *          that the hits are sorted within a layer end.
 *
 *  @param  tkr  The GLEAM representation of the hit strips + TOTs
 *
**/
void EbfTkrData::fill (const Event::TkrDigiCol &tkr, TriRowBitsTds::TriRowBits *rowbits)
{
   /*
    | Within this routine layer ends are defined as:
    |    X:  0 - 35 (even = lo strip addresses, odd = hi strips)  
    |    Y: 36 - 71 (even = lo strip addresses, odd = hi strips)
    |
    | This is a completely arbitrary definition which need only
    | be used consistently within this class.
   */

   /* Initialize the data structure for this event */
   initialize ();


   /* Protect against an empty tracker record, pretty unlikely.. */ 
   if (tkr.size ())
   {

       for (Event::TkrDigiCol::const_iterator it = tkr.begin();
            it != tkr.end();
            it++)
       {

           const Event::TkrDigi&   digi = **it;
           idents::GlastAxis::axis axis = digi.getView();
           idents::TowerId      towerId = digi.getTower();
           int                  bilayer = digi.getBilayer();
           int                lastStrip = digi.getLastController0Strip();
           unsigned int         numHits = digi.getNumHits();

//           printf("Ebf:Tower %i Layer Hit %i  View %i numHits %i ToT %i %i\n",digi.getTower().id(),digi.getBilayer(),digi.getView(),numHits,digi.getToT(0),digi.getToT(1));
           /* 
            | The layers are numbered 
            |    X layers:  0-17
            |    Y layers: 18-36
           */
           int  xy     = (axis == idents::GlastAxis::X) ? 0 : 1;
           int  ilayer = 18 * axis + bilayer;
           EbfTkrTowerData *tower = m_towers + towerId.id();


           tower->m_tots[2*ilayer]   = digi.getToT (0);
           tower->m_tots[2*ilayer+1] = digi.getToT (1);
          
           
           bool reported=false;
           /* Loop through all the hits on this bi_layer */
           for (unsigned int iHit = 0; iHit < numHits; iHit++)
           {
               int stripId = digi.getHit(iHit);

               /*
                | Check if this strip is before or after the LO HI splitpoint
               */
               EbfTkrTowerData::eLOHI loHi = (stripId > lastStrip) 
                                           ? EbfTkrTowerData::HI
                                           : EbfTkrTowerData::LO;

               /* 
                | Translate the layer numbering to a layer end numbering
                |  X layers:  0-35 
                |  Y layers: 36-71
                |  Lo end  : Even layers
                |  Hi end  : Odd  layers
               */
               int ilayerEnd = 2*ilayer + loHi;
               int     ihits =  tower->m_nhits[ilayerEnd];
               
               if(ihits<EbfTkrTowerData::MaxStripsPerLayerEnd){
                 tower->m_nhits[ilayerEnd]       = ihits + 1;
                 tower->m_data[ilayerEnd][ihits] = stripId;
                 tower->m_maps[xy][loHi]        |= (1 << bilayer);
//                 printf("Adding Hit in Tower %i and Layer %i and view %i map %5.5x\n",towerId.id(),bilayer,xy,tower->m_maps[xy][loHi]);
               } else{
			        printf("EbfTrkData:: MaxStripHits exceed in Tower %i BiLayer %i xy  %i Split %i iLayerEnd %i Hit %i StripId 0x%3.3x\n",
                         towerId.id(),bilayer,xy,loHi,ilayerEnd,iHit,stripId); 
			      }
           }
       }

      
       m_trigger = completeTkrTrg (m_towers,  
                                   sizeof (m_towers) / sizeof (*m_towers));


/*         for(int twr=0; twr<16; twr++) {
           unsigned int rwbits = rowbits->getDigiTriRowBits(twr);
           unsigned int rwbitsTrg = rowbits->getTrgReqTriRowBits(twr);
           printf("Tower %i  DigiTkr %4.4x  TrgReg %4.4x\n",twr,rwbits,rwbitsTrg);
         }
*/      
   }
   
   return;
}



/**
 *
 *  @fn     EbfTkrData::fillEncode ()
 *  @brief  Fills the hits by tower and layer end. It is assumed 
 *          that the hits are sorted within a layer end.
 *
 *  @param  tkr  The GLEAM representation of the hit strips + TOTs
 *
**/
void EbfTkrData::fillEncode (int encodeFlag, int event)
{
   /*
    | Within this routine layer ends are defined as:
    |    X:  0 - 35 (even = lo strip addresses, odd = hi strips)  
    |    Y: 36 - 71 (even = lo strip addresses, odd = hi strips)
    |
    | This is a completely arbitrary definition which need only
    | be used consistently within this class.
   */

   /* Initialize the data structure for this event */
   initialize ();
   if(encodeFlag==0) return;

   /* Loop Over towers */ 
   for (int twr=0; twr<16; twr++){
      for(int bilayer=0; bilayer<18; bilayer++){
        for(int xy=0; xy<2; xy++){
	  for(int split=0; split<2; split++){

// Set the split point and the number of hits	
           int lastStrip = 768;
           unsigned int numHits = 64;

// Setup Conditions for Small Encode (encodeFlag==2)
	   if(encodeFlag==2) numHits=1;
	   if( (bilayer != (event%18) | xy != ((event/18)%2)  | split != ((event/36)%2) ) && encodeFlag==2 ) continue;
//	   printf("Event %i bilayer %i xy %i  split %i\n",event,bilayer,xy,split);
	   
           /* 
            | The layers are numbered 
            |    X layers:  0-17
            |    Y layers: 18-36
	    /
	    /  "ilayer" must run 0-35
	    /  "bilayer must run 0-17
	    /  "xy" must be 0 or 1
           */
           int  ilayer = 18 * xy + bilayer;
           EbfTkrTowerData *tower = m_towers + twr;
         
// Don't know how to fill these at the moment. Is one of these for
// the high end and the other the low end?
           tower->m_tots[2*ilayer]   = 0;
           tower->m_tots[2*ilayer+1] = 0;
      
           bool reported=false;
           /* Loop through all the hits on this bi_layer */
           for (unsigned int iHit = 0; iHit < numHits; iHit++)
           {
	       /*
	        / Define the stripId to be encoded as follows
		/    bit   10: HI=1/LO=0
		/    bit    9:  Always Zero (Avoids going over split point)
		/    bit    8:  Tower # Even = 0 / Tower # Odd = 1
		/    bit    7:  (bi)layer # Even = 0 / Layer # Odd = 1
		/    bit    6:  x layer=0 / y layer=1 
		/    bits 0-5: hit Number  
	       */
               int stripId = iHit;
	       stripId |= xy << 6;
	       if(bilayer%2)  stripId |= 1 << 7;
	       if(twr%2)      stripId |= 1 << 8;
	       if(split==1)   stripId |= 1 << 10; 

               /*
                | Check if this strip is before or after the LO HI splitpoint
               */
               EbfTkrTowerData::eLOHI loHi = (stripId > lastStrip) 
                                           ? EbfTkrTowerData::HI
                                           : EbfTkrTowerData::LO;

               /* 
                | Translate the layer numbering to a layer end numbering
                |  X layers:  0-35 
                |  Y layers: 36-71
                |  Lo end  : Even layers
                |  Hi end  : Odd  layers
               */
               int ilayerEnd = 2*ilayer + loHi;
               int     ihits =  tower->m_nhits[ilayerEnd];
//	       printf("Tower %i BiLayer %i xy  %i Split %i iLayerEnd %i Hit %i ihit %i StripId 0x%3.3x\n",twr,bilayer,xy,split,ilayerEnd,iHit,(ihits+1),stripId);
               
               if(ihits<EbfTkrTowerData::MaxStripsPerLayerEnd){
                 tower->m_nhits[ilayerEnd]       = ihits + 1;
                 tower->m_data[ilayerEnd][ihits] = stripId;
                 tower->m_maps[xy][loHi]        |= (1 << bilayer);
               } else{
			        printf("EbfTrkData:: MaxStripHits exceed in Tower %i BiLayer %i xy  %i Split %i iLayerEnd %i Hit %i ihit %i StripId 0x%3.3x\n",twr,bilayer,xy,split,ilayerEnd,iHit,(ihits+1),stripId); 
	            }
             } // loop over hits
	   } // loop over split 
         } // loop over xy
       } // loop over bilayers
     } // loop over towers
      
       m_trigger = completeTkrTrg (m_towers,  
                                   sizeof (m_towers) / sizeof (*m_towers));

   
   return;
}


/**
 *
 *  @fn         int EbfTkrData::format (unsigned int *dst)
 *  @brief      Converts the tracker data into Event Builder Format.
 *
 *  @param  dst The destination array.
 *  @return     Pointer to the next available word in the output array.
 *
 *  It is the responsibility of the caller to ensure that enough room
 *  for a maximum record exists in the output array.  
 *
 *  Since it outputs all the towers as a contigious block, this is 
 *  just a debugging routine. The real formatter must interleave the
 *  CAL and TKR data on a tower-by=tower basis.
 *
 */
unsigned int *EbfTkrData::format (unsigned int *dst) const
{
   unsigned int itower;

   /* Loop over each tower */   
   for (itower = 0;
        itower < sizeof (m_towers) / sizeof (*m_towers);
        itower++)
   {
         /* ...and format on a tower-by-tower basis */
         dst = format (dst, itower);
   }

   return dst;
}





/**
 *
 *  @fn     unsigned int *EbfTkrData::format (unsigned int *dst,
 *                                           int         itower) const
 *  @brief  Formats the TKR data for one tower into Event Builder Format
 *
 *  @param    dst The destionation array.
 *  @param itower The tower to format.
 *  @return       Pointer to the next available word in the destination array.
 *
 */
unsigned int *EbfTkrData::format (unsigned int *dst, int itower) const
{
   /*
    |
    | The numbering scheme is such that
    |
    |   X_0 lo =  0       Y_0 lo = 36
    |   X_0 hi =  1       Y_0 hi = 37
    |
    |   X_1 lo =  2       Y_1 lo = 38
    |   X_1 hi =  3       Y_1 hi = 39
    |
    |   X_2 lo =  4       Y_2 lo = 40
    |   X_2 hi =  5       Y_2 hi = 41
    |
    |   X_3 lo =  6       Y_3 lo = 42
    |   X_3 hi =  7       Y_3 hi = 43
    |
    |   X_4 lo =  8       Y_4 lo = 44
    |   X_4 hi =  9       Y_4 hi = 45
    |
    |   X_5 lo = 10       Y_5 lo = 36
    |   X_5 hi = 11       Y_5 hi = 37
    |
    |   X_6 lo = 12       Y_6 lo = 38
    |   X_6 hi = 13       Y_6 hi = 39
    |
    |   X_7 lo = 14       Y_7 lo = 40
    |   X_7 hi = 15       Y_7 hi = 41
    |
    |   X_8 lo = 16       Y_8 lo = 42
    |   X_8 hi = 17       Y_8 hi = 43
    |
    |   X_9 lo = 18       Y_9 lo = 44
    |   X_9 hi = 19       Y_9 hi = 45
    |
    |   X10 lo = 20       Y10 lo = 46
    |   X10 hi = 21       Y10 hi = 47
    |
    |   X11 lo = 22       Y11 lo = 48
    |   X11 hi = 23       Y11 hi = 49
    |
    |   X12 lo = 24       Y12 lo = 50
    |   X12 hi = 25       Y12 hi = 51
    |
    |   X13 lo = 26       Y13 lo = 52
    |   X13 hi = 27       Y13 hi = 53
    |
    |   X14 lo = 28       Y14 lo = 64
    |   X14 hi = 29       Y14 hi = 65
    |
    |   X15 lo = 30       Y15 lo = 66
    |   X15 hi = 31       Y15 hi = 67
    |
    |   X16 lo = 32       Y16 lo = 68
    |   X16 hi = 33       Y16 hi = 69
    |
    |   X17 lo = 34       Y17 lo = 70
    |   X17 hi = 35       Y17 hi = 71
   */

   /* 
    | Using the above information, fill a data structure giving the number
    | of the first layer on a controller
   */

   static const int Begin[EbfTkrTowerData::NumCables] = { 2,  0,  3,  1,
                                                         38, 36, 39, 37  };
   

   unsigned int  accepts[EbfTkrTowerData::NumCables];
   int                                    nlayerEnds;
   int                                          boff;
   const EbfTkrTowerData   *tower = m_towers + itower;   
   int           tots[EbfTkrTowerData::NumLayerEnds];
       
   /* 
    | First compute the accept bit list for the cable ends 
    | The output cable order is
    |      6. X odd  lo ( 2,  6, 10, 14, 18, 22, 26, 30, 34)
    |      3. X even lo ( 0,  4,  8, 12, 16, 20, 24, 28, 32)
    |      7. X odd  hi ( 3,  7, 11, 15, 19, 23, 27, 31, 35)
    |      2. X even hi ( 1,  5,  9, 13, 17, 21, 25, 29, 33)
    |      5. Y even hi (38, 42, 46, 50, 54, 58, 62, 66, 70)
    |      0. Y even lo (36, 40, 44, 48, 52, 56, 60, 64, 68)
    |      4. Y odd  hi (39, 43, 47, 51, 55, 59, 63, 67, 71)
    |      1. Y odd  lo (37, 41, 45, 49, 53, 57, 61, 65, 69)
   */
   nlayerEnds = 0;
   for (unsigned int ocable = 0; ocable<EbfTkrTowerData::NumCables; ocable++)
   {
       /*
        | mlayer is the accept mask for the current layer. The lowest layer
        | is placed in the MSB of the layer accept mask.
       */
       int mlayer       = (1 << (EbfTkrTowerData::NumLayerEndsPerCable - 1));
       int ilayerEndMin = Begin[ocable];
       int ilayerEndMax = ilayerEndMin 
                        + 4*EbfTkrTowerData::NumLayerEndsPerCable;
       int         list = 0;

       /* Pack the 9 layer accept bits for this cable */
       for (int ilayerEnd = ilayerEndMin;
            ilayerEnd     < ilayerEndMax;
            ilayerEnd    += 4, mlayer >>= 1)
       {
           /* 
            | The number of layer ends with hits must be kept track of
            | so that one can pack one TOT per hit layer end. The tots
            | are sorted into an array so they can by jammed out in
            | the correct order.
           */
           if (tower->m_nhits[ilayerEnd]) 
           {
               tots[nlayerEnds] = tower->m_tots[ilayerEnd];
               nlayerEnds      += 1;
               list            |= mlayer;
           }
       }

       accepts[ocable] = list;
   }

   boff = packAccepts (dst, accepts);

   int totalHits[EbfTkrTowerData::NumCables];
   
   /* Now pack up the hits on each cable */
   for (unsigned int ocable = 0; ocable<EbfTkrTowerData::NumCables; ocable++)
   {
       totalHits[ocable] = 0;
       int ilayerEndMin = Begin[ocable];
       int ilayerEndMax = ilayerEndMin 
                        + 4*EbfTkrTowerData::NumLayerEndsPerCable;
       
//       printf("Tower %2i Cable %2i",itower,ocable);
       /* Pack the layer hits on the cable */
       for (int ilayerEnd = ilayerEndMin;
            ilayerEnd     < ilayerEndMax;
            ilayerEnd    += 4)
       {
           boff = packLayerHits (dst, 
                                 boff, 
                                 tower->m_data[ilayerEnd],
                                 tower->m_nhits[ilayerEnd],
                                 totalHits[ocable]);
           totalHits[ocable] +=tower->m_nhits[ilayerEnd]; 
//           printf(" %4i",tower->m_nhits[ilayerEnd]);
       }
//       if(totalHits[ocable]>100) printf("EbfTrkData:: Tower %i Cable %i Hits = %4i \n",itower,ocable,totalHits[ocable]); 
       if(totalHits[ocable]>EbfTkrTowerData::MaxStripsPerCable) {
            printf("EbfTrkData:: Tower %i Cable %i has too many hits %i \n",itower,ocable,totalHits[ocable]);
//            tower->m_Truncated[ocable] = true;
       } 
   }

   boff = packTots (dst, boff, nlayerEnds, tots);
   

   /* 
    | Return the address of the next available 32 bit word.
    | If the current bit offset stops within a 32 bit word, then round
    | to the next word, else if the bit offset is an exact multiple of
    | 32, then the no need to round, since no bits are used from this word.
   */
   return &dst[boff >> 5] + ((boff & 0x1f) != 0);

}


void EbfTkrData::parseInput(unsigned int *contrib, unsigned int itower, unsigned int lcbWords)
{

   static const int Begin[EbfTkrTowerData::NumCables] = { 2,  0,  3,  1,
                                                         38, 36, 39, 37  };

   unsigned int *data = contrib;
   bool debug = false;
   
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
   
   if(itower==9) initialize();
   
   
// Get the tower
   EbfTkrTowerData *tower = m_towers + itower;
   
   
// Peel off header words   
   unsigned int LCB_Header0 = *data++;
   unsigned int LCB_Header1 = *data++;

// Parse out header words   
   unsigned int startbit  = (LCB_Header1 >> 31) & 0x1;
              m_calStrobe = (LCB_Header1 >> 30) & 0x1;
   unsigned int sequence0 = (LCB_Header1 >> 28) & 0x3;
   unsigned int trigAck   = (LCB_Header1 >> 27) & 0x1;
   unsigned int range4    = (LCB_Header1 >> 26) & 0x1;
   unsigned int zeroSup   = (LCB_Header1 >> 25) & 0x1; 
   unsigned int marker    = (LCB_Header1 >> 22) & 0x7;
   unsigned int ErrCont   = (LCB_Header1 >> 21) & 0x1;
   unsigned int DiagCont  = (LCB_Header1 >> 20) & 0x1;
   unsigned int sequence1 = (LCB_Header1 >> 1) & 0x7fff;
   unsigned int ParErr    = (LCB_Header1 & 0x1);

    if(debug) printf("LCB_Header1 0x%8.8x\n",LCB_Header1);
    if(debug) printf("startbit %1i  calStrobe %1i seq %1i trigAck %1i 4Rng %1i ZSup %1i \nMkr %1i ErrCnt %1i DiaCnt %1i ParEr %1i Seq %i\n",
            startbit,m_calStrobe,sequence0,trigAck,range4,zeroSup,marker,ErrCont,DiagCont,ParErr,sequence1);
                
// Get the number of Calorimeter hits and skip over them
   unsigned int calHits_Wd = *data++;
   unsigned int calHits = 0;
   for(int i=0; i<8; i++) calHits += (calHits_Wd >> i*4) & 0xf;   
   if(debug) printf("Hits in Cal Layers: %i\n",calHits);

// Advance past the calorimeter information
   int nmWds = range4 == 1 ? calHits*4 : calHits;
   for(int i=0; i<nmWds; i++) data++;

// Grab the cable Accept bits
   unsigned int cableWd0 = *data++;
   unsigned int cableWd1 = *data++;
   unsigned int cableWd2 = *data;  //Don't advance pointer because we have strip info on this word

// Now Parse Accept bits:
   unsigned int Accept[8];
   Accept[0] = (cableWd0 >>23) & 0x1ff;   
   Accept[1] = (cableWd0 >>14) & 0x1ff;
   Accept[2] = (cableWd0 >>5) & 0x1ff;
   Accept[3] = ((cableWd0 & 0x1f) << 4) | ((cableWd1 >>28) & 0xf) ;
   Accept[4] = (cableWd1 >>19) & 0x1ff;
   Accept[5] = (cableWd1 >>10) & 0x1ff;
   Accept[6] = (cableWd1 >>1) & 0x1ff;
   Accept[7] = ((cableWd1 & 0x1) << 8) | ((cableWd2 >>24) & 0xff) ;
   if(debug) for(int i=0; i<8; i++) printf("Accept[%i]: 0x%8.8x\n",i,Accept[i]);

// Now Start grabbing the hits
   int bitOff = 12; // Starting point after cable accepts.
   unsigned int stripID = 0;
   unsigned int nTOTs =0;
   unsigned int TOTLayer[EbfTkrTowerData::NumLayerEnds];
   
   for(int cable=0; cable<8; cable++) {

// Set the starting layerEnd for this cable
      unsigned int ilayerEnd = Begin[cable];

      unsigned int accList = Accept[cable];
      while(accList != 0) {

// first layer is in most sig. bit of accept list peel off starting at
// the top bits      
       if( (accList & 0x100) != 0) {
         TOTLayer[nTOTs] = ilayerEnd;         
         nTOTs++;
         if(debug) printf("Cable %i accList 0x%8.8x\n",cable,accList);
         bool LastHit = false; 
         while(!LastHit) {
           
// Figure out where we need to go for strip ID
           if(debug) printf("BitOff = %i\n",bitOff);
           if(bitOff >= 0) {
              stripID = (*data >> bitOff) & 0xfff;
              if(debug) printf("NonNeg: bitOff %i  stripID 0x%8.8x\n",bitOff,stripID);
           } else {
              stripID =  (*data++ & (0xfff >> -bitOff)) << -bitOff;
              if(debug) printf("Neg1: bitOff %i  stripID 0x%8.8x\n",bitOff,stripID);
              bitOff += 32;
              stripID |= (*data >> bitOff); 
              if(debug) printf("Neg2: bitOff %i  stripID 0x%8.8x\n",bitOff,stripID);
           }
           
// Set strip ID for this layer
           int ihits = tower->m_nhits[ilayerEnd];
           tower->m_nhits[ilayerEnd] = ihits+1;
           tower->m_data[ilayerEnd][ihits] = stripID;
                      
// shift the bit offset                      
           bitOff -= 12;

// Check to see if this is the last strip for this layer                         
           if( (stripID&0x800) !=0 ) LastHit = true;
           
         } //while more hits on this layer
        } //if this layer has an accept bit set
       accList <<= 0x1;

// Set layer end for this accept list
       ilayerEnd += 4;               
       
      }//while accept bits
   } // loop over cables   
        
// Now Get the TOT
   unsigned int TOT = 0;
   if(debug) printf("Number of TOTs %i\n",nTOTs);
   for(int i=0; i<nTOTs; i++) {
      if(debug) printf("start bitOff = %i", bitOff);
      if(bitOff==-4) bitOff = 0; //First TOT is at end of 32 bit word
      if(bitOff==-12) { bitOff = 24; data++; } //32 bit word is finished by padding..advance to next word
      if(bitOff%8 != 0) bitOff -=4; //middle of 32 bits word...get to byte boundary
      if(debug) printf("End bitOff = %i\n",bitOff);
      
      if(bitOff >= 0) {
         TOT = (*data >> bitOff) & 0xff;
         if(debug) printf("NonNeg: bitOff %i  TOT 0x%8.8x\n",bitOff,TOT);
      } else {
         TOT =  (*data++ & (0xff >> -bitOff)) << -bitOff;
         if(debug) printf("Neg1: bitOff %i  TOT 0x%8.8x\n",bitOff,TOT);
         bitOff += 32;
         TOT |= (*data >> bitOff); 
         if(debug) printf("Neg2: bitOff %i  TOT 0x%8.8x\n",bitOff,TOT);
      }

// Set the TOT. Ebf Has dropped two least sig. bits...must scale up but
// resolution is lost compared to original mc digi's
      tower->m_tots[TOTLayer[i]]= TOT << 2;

// Move the bit offset            
      bitOff -= 8;          
   }

   return;
}
  
/**
 *
 *  @fn     void EbfTkrData::print ()
 *  @brief  Diagnostic print routine 
 * 
 *  Routine to print the hits on each layer end of each tower. This
 *  is mainly used to debug the filling routine
 *
**/
void EbfTkrData::print()
{
   EbfTkrTowerData   *tower;
   unsigned int      itower;

   std::printf ("TKR: Trigger = %8.8x\n", m_trigger);

   /* Loop over each tower */
   for (itower = 0, tower = m_towers;
        itower < sizeof (m_towers) / sizeof (*m_towers); 
        itower++, tower++)
   {
       unsigned int ilayerEnd;

       std::printf ("\nTwr:%1.1x X:(%8.8x | %8.8x) Y:(%8.8x | %8.8x)\n", 
                    itower, 
                    tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::LO],
                    tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::HI],
                    tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::LO],
                    tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::HI]);
       

       /* Loop over each layer end of the current tower */ 
       for (ilayerEnd = 0; 
            ilayerEnd < sizeof (tower->m_nhits) / sizeof (*tower->m_nhits); 
            ilayerEnd++)
       {
           int                    ihits;
           int                     ncol;
           int                  nmargin;
           int                   ilayer;
           int                    nhits;
           unsigned short int *stripIds;


           /* 
            | Check if there any hits on this layer.
            | If none, then on to the next layer.
           */
           nhits = tower->m_nhits[ilayerEnd];
           if (nhits == 0) continue;

           /*
            | Have some hits.
            |  1. Compute a pointer to the strip ids on this layer.
            |  2. Change the layer end into a just a layer number
            |  3. Print the preamble for the hits
           */
           stripIds = tower->m_data[ilayerEnd];       
           ilayer   = ilayerEnd < 36 ? ilayerEnd : ilayerEnd - 36;
           nmargin  = std::printf ("L%2d%c%c", 
                                   ilayer/2,
                                   (ilayer&1)  ? 'h' : 'l',
                                   ilayerEnd < 36 ? 'X' : 'Y');
           ncol     = nmargin;
           
       

           /* Loop over all the strip hits on this layer end */
           for (ihits = 0; ihits < nhits; ihits++)
           {
               unsigned short int stripId = stripIds[ihits];
           
               if (ncol > 75) ncol = std::printf ("\n%*c", nmargin, ' ') - 1;
               ncol += std::printf (" %4d:%3.3x", stripId, stripId);
           }
           std::printf ("\n");
       }
   }

   std::printf ("\n");
   return;
}



/**
 *
 *  @fn unsigned int completeTkrTrg (const EbfTkrTowerData *tower,
 *                                  unsigned int         ntowers)
 *  @brief Calculates the tracker trigger for each of the towers.
 *
 *  @param   tower The tracker tower records
 *  @param ntowers The number of towers to process.
 *  @return        A bit mask indicating which towers had a tracker trigger
 *                 LSB = Tower 0, MSB = Tower 15
 *
 */
static unsigned int completeTkrTrg (const EbfTkrTowerData *tower, 
                                    unsigned int           ntowers)
{
   unsigned int    trigger;
   unsigned int     itower;

   trigger = 0;

   /* This loop forms the basic 3-in-a-row trigger for each tower */
   for (itower = 0; itower < ntowers; itower++, tower++)
   {
       unsigned int   x;
       unsigned int   y;
       unsigned int  xy;
       unsigned int tir;
       
       x = tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::LO] 
         | tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::HI];

       y = tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::LO]
         | tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::HI];

   
       xy  = x & y;
       tir = (xy >> 2) & (xy >> 1) & (xy >> 0);
       if (tir) trigger |= (1 << itower);

//               printf("Maps twr %i: xh %4.4x xl %4.4x  yh %4.4x  yl %4.4x  x=%4.4x y=%4.4x  xy=%4.4x tir=%4.4x\n",
//               itower,tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::HI],tower->m_maps[EbfTkrTowerData::X][EbfTkrTowerData::LO],
//               tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::HI],tower->m_maps[EbfTkrTowerData::Y][EbfTkrTowerData::LO],
//               x,y,xy,tir);

   }
   
   return trigger;
}






/**
 *
 * @fn int packAccepts (unsigned int                                    *dst,
 *                      const unsigned int accepts[EbfTkrTowerData::NumCables])
 * @brief       Packs the 8 cable accept lists into the 72 bit output accept
 *              list
 *
 * @param   dst The destination array to receive the packed values.
 * @param  boff The current bit offset into the destination array. Bit 0
 *               = MSB, consistent with big endian addressing.
 * @param accepts The array of the 8 cable layer accepts
 * @return        The bit offset to store the next data values at. Given
 *                that this is big endian packing, this bit offset is
 *                from the top, ie Bit 0 = MSB
 *
**/
static int 
       packAccepts (unsigned int                                      *dst,
                    const unsigned int accepts[EbfTkrTowerData::NumCables])
{
   /*
    |      3322222222221111111111
    |      10987654321098765432109876543210
    |      --------------------------------
    |      00000000011111111122222222233333
    |      33334444444445555555556666666667
    |      77777777
   */

   /* 
    | First of the accept list words, 
    | note only 5 MSB  of the 3rd accept word fits
   */  
   dst[0] = (accepts[0] << (32- 9))   /*   0 - 8    LSB = 23  # of bits = 9 */
          | (accepts[1] << (32-18))   /*   9 - 17   LSB = 14  # of bits = 9 */
          | (accepts[2] << (32-27))   /*  18 - 26   LSB =  5  # of bits = 9 */
          | (accepts[3] >> (36-32));  /*  27 - 31   LSB =  0  # of bits = 5 */


   /* 
    | Second of the accept list words,
    | note only 4 LSB of the 3rd accept word
    | and  only 1 MSB of the 7th accept fits 
   */
   dst[1] = (accepts[3] << (32-4))     /*  32 - 35  LSB = 28  # of bits = 4 */ 
          | (accepts[4] << (32-4- 9))  /*  36 - 44  LSB = 19  # of bits = 9 */ 
          | (accepts[5] << (32-4-18))  /*  45 - 53  LSB = 10  # of bits = 9 */ 
          | (accepts[6] << (32-4-27))  /*  54 - 62  LSB =  1  # of bits = 9 */
          | (accepts[7] >> (4+36-32)); /*  63 - 63  LSB =  0  # of bits = 1 */


   /* Good till the last byte */
   dst[2] = (accepts[7] << (32-8));    /*  64 - 71  LSB = 24  # of bits = 8 */

   return EbfTkrTowerData::NumLayerEnds;
}





/**
 *
 *  @fn    int packLayerHits (unsigned int                  *dst, 
 *                            int                           boff, 
 *                            const unsigned short int *stripIds,
 *                            int                          nhits)
 *
 *  @brief  Packs the hits for 1 layer into the destination array.
 *
 *  @param      dst The destination array to receive the packed values.
 *  @param     boff The current bit offset into the destination array.
 *  @param stripIds The array of strip ids to pack,
 *  @param    nhits The number of strip ids to pack.
 *  @return         The next bit offset into the destination array.
 *
**/
static int packLayerHits (unsigned int                  *dst, 
                          int                           boff, 
                          const unsigned short int *stripIds,
                          int                          nhits,
                          int                       cableHits)
{
   /*
    |  Loop over the hits, packing the strip ids 12  bits at a time
   */
   
   // Keep track of cable hits for truncation of Ebf 
   int totcableHits = cableHits;
   
//   printf("Packing the strip IDs: nhits = %i\n",nhits);
   while (nhits > 0 && (totcableHits++)<EbfTkrTowerData::MaxStripsPerCable)
   {
       int              idx;
       int              bdx;
       unsigned int stripId;
               
       stripId = *stripIds++;

       /* If at the end of a layer, add the terminator bit */
       if (--nhits <= 0 || (totcableHits+1)==EbfTkrTowerData::MaxStripsPerCable) stripId |= 0x800;


       /* 
        | Determine where the current offset within a 32 bit word
        | boff = absolute bit position of LSB of the strip of the 
        |        previous strip with MSB bit0 (after addition of 12).
        | idx  = index of the 32-bit word holding the low address
        | bdx  = right shift necessary to position the strip address
        |
        |                33222222 22221111 111111 
        |                10987654 32109876 54321098 76543210
        |  boff  bdx    +--------+--------+--------+--------+
        |     0   20     ssssssss ssss
        |    12    8                  ssss ssssssss
        |    24   -4                                hhhhhhhh llll
        |    36   16     hhhhssss ssssssss
        |    48    4                       ssssssss ssss
        |    60   24                                    hhhh llllllll
        |    72   12     llllllll ssssssss ssss
        |    84    0                           ssss ssssssss 
        |    96   20     ssssssss ssss
        |
       */
       idx   = boff >> 5;
       bdx   = ((boff ^ 0x1f) & 0x1f) - 11;
       boff += 12;
       

       if (bdx == 32-12) 
       {
           /* At the beginning, new word, just store the strip id */ 
           dst[idx] = stripId << 20;
       }
       else if (bdx >= 0)
       {
           /* Strip fits, add to the previous word */
           dst[idx] |= stripId << bdx;
       }
       else 
       {
           /* Split the strip address across the boundary              */
           dst[idx]   |= stripId >>   (-bdx);      /* Get the hi piece */
           dst[idx+1]  = stripId << (32+bdx);      /* Get the lo piece */
       }
   }

   return boff;
}





/**
 *
 *  @fn   int packTots (unsigned int *dst,
 *                     int           boff,
 *                     int          ntots)
 *
 *  @brief Packs the tracker TOT values into the destination array.
 *
 *  @param   dst The destination array to receive the packed values.
 *  @param  boff The current bit offset into the destination array.
 *  @param ntots The number of TOTs to pack
 *  @return      The next bit offset
 *
 */
static int packTots (unsigned int                  *dst,
                     int                           boff,
                     int                          ntots,
                     const int                    *tots)
{
   int itot;
	
   /*
    | Now pack the TOTs 
    | First round the bit offset to the nearest byte
   */
//   printf("Packing TOTs: ntots = %i initial boff %i\n",ntots,boff);
   boff = (boff + 7) & ~0x7;
//   printf("Revised boff %i\n",boff);
   for (itot = 0; itot < ntots; itot++)
   {
       int bdx;
       int tot;
       
       tot   = *tots++;

// The TEM truncates the TOT value to 8 bits by dropping the
// least significant 2 bits.
       tot  = (tot>>2) & 0xff;
       
       bdx   = ((boff ^ 0x1f) & 0x1f) - 7;
        
//       printf("Storing TOT %i boff %i  bdx %i  tot 0x%8.8x \n",itot,boff,bdx,tot);
       if   (bdx == 24) dst[boff >> 5]  = tot <<  24;
       else             dst[boff >> 5] |= tot << bdx;
       
       boff += 8;
   }

   return boff;
}




#define ConfigCompare_cxx
#include "configData/ConfigCompare.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>

void ConfigCompare::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ConfigCompare.C
//      Root > ConfigCompare t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

Bool_t 
ConfigCompare::compare(ConfigCompare& other, Long64_t entry, Bool_t fullCompare, Bool_t onlyOne) {
  
  LoadTree(entry); GetEntry(entry);
  other.LoadTree(entry); other.GetEntry(entry);

  Bool_t latch = kTRUE;

  // first singletons
  latch &= compareSingleton("aem_configuration",aem_configuration,other.aem_configuration);
  latch &= compareSingleton("trgseq",trgseq,other.trgseq);
  latch &= compareSingleton("window_width",window_width,other.window_width);
  latch &= compareSingleton("window_open_mask",window_open_mask,other.window_open_mask);
  latch &= compareSingleton("periodic_prescale",periodic_prescale,other.periodic_prescale);
  latch &= compareSingleton("periodic_limit",periodic_limit,other.periodic_limit);
  latch &= compareSingleton("periodic_free",periodic_free,other.periodic_free);
  latch &= compareSingleton("periodic_onepps",periodic_onepps,other.periodic_onepps);
  latch &= compareSingleton("acd_cno",acd_cno,other.acd_cno);
  latch &= compareSingleton("tower_busy",tower_busy,other.tower_busy);
  latch &= compareSingleton("external",external,other.external);
  latch &= compareSingleton("configuration",configuration,other.configuration);

  // simple vectors
  //compareV1("hold_delay",hold_delay,other.hold_delay,12);
  latch &= compareV1("hitmap_delay",hitmap_delay,other.hitmap_delay,12,onlyOne);  
  latch &= compareV1("hitmap_width",hitmap_width,other.hitmap_width,12,onlyOne);
  latch &= compareV1("veto_delay",veto_delay,other.veto_delay,12,onlyOne);
  latch &= compareV1("veto_width",veto_width,other.veto_width,12,onlyOne);

  latch &= compareV1("data_masks",data_masks,other.data_masks,16,onlyOne);
  latch &= compareV1("tkr_trgseq",tkr_trgseq,other.tkr_trgseq,16,onlyOne);
  latch &= compareV1("cal_trgseq",cal_trgseq,other.cal_trgseq,16,onlyOne);
 
  latch &= compareV1("conditions",conditions,other.conditions,32,onlyOne);
  latch &= compareV1("engine",engine,other.engine,16,onlyOne);
  latch &= compareV1("towers",towers,other.towers,4,onlyOne);
  latch &= compareV1("tiles",tiles,other.tiles,12,onlyOne);
  latch &= compareV1("r",r,other.r,54,onlyOne);

  // arrays
  latch &= compareV2("tci_dac",(UShort_t*)tci_dac,(UShort_t*)other.tci_dac,12,18,onlyOne);

  if ( fullCompare ) {
    latch &= compareV2("pha_threshold",(UShort_t*)pha_threshold,(UShort_t*)other.pha_threshold,12,18,onlyOne);
    latch &= compareV2("bias_dac",(UShort_t*)bias_dac,(UShort_t*)other.bias_dac,12,18,onlyOne);
    latch &= compareV2("hld_dac",(UShort_t*)hld_dac,(UShort_t*)other.hld_dac,12,18,onlyOne);
    latch &= compareV2("veto_dac",(UShort_t*)veto_dac,(UShort_t*)other.veto_dac,12,18,onlyOne);   
    latch &= compareV2("veto_vernier",(UShort_t*)veto_vernier,(UShort_t*)other.veto_vernier,12,18,onlyOne);
  }

  latch &= compareV2("low",(UInt_t*)low,(UInt_t*)other.low,16,36,onlyOne);   
  latch &= compareV2("high",(UInt_t*)high,(UInt_t*)other.high,16,36,onlyOne); 
  
  latch &= compareV2("tcc_trg_align",(UInt_t*)tcc_trg_align,(UInt_t*)other.tcc_trg_align,16,8,onlyOne);   
  latch &= compareV2("tcc_configuration",(UInt_t*)tcc_configuration,(UInt_t*)other.tcc_configuration,16,8,onlyOne);   
  latch &= compareV2("input_mask",(UInt_t*)input_mask,(UInt_t*)other.input_mask,16,8,onlyOne);  
 
  latch &= compareV2("ccc_trg_alignment",(UInt_t*)ccc_trg_alignment,(UInt_t*)other.ccc_trg_alignment,16,4,onlyOne);   
  latch &= compareV2("ccc_configuration",(UInt_t*)ccc_configuration,(UInt_t*)other.ccc_configuration,16,4,onlyOne);   
  latch &= compareV2("layer_mask_0",(UInt_t*)layer_mask_0,(UInt_t*)other.layer_mask_0,16,4,onlyOne);   
  latch &= compareV2("layer_mask_1",(UInt_t*)layer_mask_1,(UInt_t*)other.layer_mask_1,16,4,onlyOne);   

  latch &= compareV2("crc_dac",(UInt_t*)crc_dac,(UInt_t*)other.crc_dac,16,4,onlyOne);   
  latch &= compareV2("delay_1",(UInt_t*)delay_1,(UInt_t*)other.delay_1,16,4,onlyOne);   
  latch &= compareV2("delay_2",(UInt_t*)delay_2,(UInt_t*)other.delay_2,16,4,onlyOne);   
  latch &= compareV2("delay_3",(UInt_t*)delay_3,(UInt_t*)other.delay_3,16,4,onlyOne);   
  
  latch &= compareV3("trig_enable",(ULong64_t*)trig_enable,(ULong64_t*)other.trig_enable,16,36,24,onlyOne);
  latch &= compareV3("injection",(UInt_t*)injection,(UInt_t*)other.injection,16,36,24,onlyOne);
  
  if ( fullCompare ) {  
    latch &= compareV3("data_mask",(ULong64_t*)data_mask,(ULong64_t*)other.data_mask,16,36,24,onlyOne);
    latch &= compareV3("calib_mask",(ULong64_t*)calib_mask,(ULong64_t*)other.calib_mask,16,36,24,onlyOne);
    latch &= compareV3("threshold",(UInt_t*)threshold,(UInt_t*)other.threshold,16,36,24,onlyOne);
  }

  latch &= compareV3("trc_csr",(UInt_t*)trc_csr,(UInt_t*)other.trc_csr,16,8,9,onlyOne);
  
  latch &= compareV4("config_0",(UInt_t*)config_0,(UInt_t*)other.config_0,16,4,4,12,onlyOne);
  latch &= compareV4("config_1",(UInt_t*)config_1,(UInt_t*)other.config_1,16,4,4,12,onlyOne);
  latch &= compareV4("ref_dac",(UInt_t*)ref_dac,(UInt_t*)other.ref_dac,16,4,4,12,onlyOne);
  
  if ( fullCompare ) {  
    latch &= compareV4("log_acpt",(UInt_t*)log_acpt,(UInt_t*)other.log_acpt,16,4,4,12,onlyOne);
    latch &= compareV4("fle_dac",(UInt_t*)fle_dac,(UInt_t*)other.fle_dac,16,4,4,12,onlyOne);
    latch &= compareV4("fhe_dac",(UInt_t*)fhe_dac,(UInt_t*)other.fhe_dac,16,4,4,12,onlyOne);
    latch &= compareV4("rng_uld_dac",(UInt_t*)rng_uld_dac,(UInt_t*)other.rng_uld_dac,16,4,4,12,onlyOne);
  }
  return latch;
} 
  
Bool_t ConfigCompare::compareSingleton(const char* name, ULong64_t val1, ULong64_t val2) {
  //static const ULong64_t intMask(0x00000000FFFFFFFF);
  if ( val1 != val2 ) {
    printf("%s ",name);
    std::cout << std::hex << val1 << ' ' << val2 << std::dec << std::endl;
    return kFALSE;
  }
  return kTRUE;
}
Bool_t ConfigCompare::compareSingleton(const char* name, UInt_t val1, UInt_t val2) {
  if ( val1 != val2 ) {
    printf("%s 0x%0x 0x%0x\n",name,val1,val2);
    return kFALSE;
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareSingleton(const char* name,  Char_t val1, Char_t val2) {
  if ( val1 != val2 ) {
    printf("%s 0x%0x 0x%0x\n",name,val1,val2);
    return kFALSE;
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV1(const char* name, UInt_t* val1, UInt_t* val2, int size, Bool_t onlyOne) {
  for (int i(0); i < size; i++ ) {
    if ( val1[i] != val2[i] ) {
      printf("%s[%i] 0x%0x 0x%0x\n",name,i,val1[i],val2[i]);
      if (onlyOne) return kFALSE;
    }
  }
  return kTRUE;
};

Bool_t ConfigCompare::compareV1(const char* name, UShort_t* val1, UShort_t* val2, int size, Bool_t onlyOne) {
  for (int i(0); i < size; i++ ) {
    if ( val1[i] != val2[i] ) {
      printf("%s[%i] 0x%0x 0x%0x\n",name,i,val1[i],val2[i]);
     if (onlyOne)  return kFALSE;
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV2(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, Bool_t onlyOne) {
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      int idx = i*size2 + j;
      if ( val1[idx] != val2[idx] ) {
	printf("%s[%i,%i] 0x%0x 0x%0x\n",name,i,j,val1[idx],val2[idx]);
	if (onlyOne) return kFALSE;
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV2(const char* name, UShort_t* val1, UShort_t* val2, int size1, int size2, Bool_t onlyOne){
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      int idx = i*size2 + j;
      if ( val1[idx] != val2[idx] ) {
	printf("%s[%i,%i] 0x%0x 0x%0x\n",name,i,j,val1[idx],val2[idx]);
	if (onlyOne) return kFALSE;
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV3(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, int size3, Bool_t onlyOne){
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {
	int idx = i*size2*size3 + j*size3 +k;
	if ( val1[idx] != val2[idx] ) {
	  printf("%s[%i,%i,%i] 0x%0x 0x%0x\n",name,i,j,k,val1[idx],val2[idx]);
	  if (onlyOne) return kFALSE;
	}
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV3(const char* name, ULong64_t* val1, ULong64_t* val2, int size1, int size2, int size3, Bool_t onlyOne){
  //static const UInt_t intMask(0x00000000FFFFFFFF);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {
	int idx = i*size2*size3 + j*size3 +k;
	if ( val1[idx] != val2[idx] ) {
	  printf("%s[%i,%i,%i] ",name,i,j,k);
	  std::cout << std::hex << val1[idx] << ' ' << val2[idx] << std::dec << std::endl;
	  if (onlyOne) return kFALSE;
	}
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV4(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, int size3, int size4, Bool_t onlyOne){
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {
	for (int l(0); l < size4; l++ ) {
	  int idx = i*size2*size3*size4 + j*size3*size4 + k*size4 + l;
	  if ( val1[idx] != val2[idx] ) {
	    printf("%s[%i,%i,%i,%i] 0x%0x 0x%0x\n",name,i,j,k,l,val1[idx],val2[idx]);
	    if (onlyOne) return kFALSE;
	  }
	}
      }
    }
  }
  return kTRUE;
}
  


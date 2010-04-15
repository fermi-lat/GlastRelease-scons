#define ConfigCompare_cxx
#include "configData/ConfigCompare.h"
#include "configData/ConfigBase.h"
#include <TH2.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH1S.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TList.h>
#include <stdio.h>
#include <iostream>


Bool_t 
ConfigCompare::compare(ConfigBase& one, ConfigBase& other, TList& l, 
		       Bool_t fullCompare, Bool_t onlyOne) {
  
  Bool_t latch = kTRUE;

  // first singletons
  latch &= compareSingleton("aem_configuration",one.aem_configuration,other.aem_configuration);
  latch &= compareSingleton("trgseq",one.trgseq,other.trgseq);
  latch &= compareSingleton("window_width",one.window_width,other.window_width);
  latch &= compareSingleton("window_open_mask",one.window_open_mask,other.window_open_mask);
  latch &= compareSingleton("periodic_prescale",one.periodic_prescale,other.periodic_prescale);
  latch &= compareSingleton("periodic_limit",one.periodic_limit,other.periodic_limit);
  latch &= compareSingleton("periodic_free",one.periodic_free,other.periodic_free);
  latch &= compareSingleton("periodic_onepps",one.periodic_onepps,other.periodic_onepps);
  latch &= compareSingleton("acd_cno",one.acd_cno,other.acd_cno);
  latch &= compareSingleton("tower_busy",one.tower_busy,other.tower_busy);
  latch &= compareSingleton("external",one.external,other.external);
  latch &= compareSingleton("configuration",one.configuration,other.configuration);

  // simple vectors
  latch &= compareV1("hold_delay",one.hold_delay,other.hold_delay,12);
  latch &= compareV1("hitmap_delay",one.hitmap_delay,other.hitmap_delay,12,onlyOne);  
  latch &= compareV1("hitmap_width",one.hitmap_width,other.hitmap_width,12,onlyOne);
  latch &= compareV1("hitmap_deadtime",one.hitmap_deadtime,other.hitmap_deadtime,12,onlyOne);
  latch &= compareV1("veto_delay",one.veto_delay,other.veto_delay,12,onlyOne);
  latch &= compareV1("veto_width",one.veto_width,other.veto_width,12,onlyOne);
  latch &= compareV1("adc_tacq",one.adc_tacq,other.adc_tacq,12,onlyOne);

  latch &= compareV1("data_masks",one.data_masks,other.data_masks,16,onlyOne);
  latch &= compareV1("tkr_trgseq",one.tkr_trgseq,other.tkr_trgseq,16,onlyOne);
  latch &= compareV1("cal_trgseq",one.cal_trgseq,other.cal_trgseq,16,onlyOne);

  latch &= compareV1("cal_in_mask",one.cal_in_mask,other.cal_in_mask,16,onlyOne);
  latch &= compareV1("tkr_layer_enable_0",one.tkr_layer_enable_0,other.tkr_layer_enable_0,16,onlyOne);
  latch &= compareV1("tkr_layer_enable_1",one.tkr_layer_enable_1,other.tkr_layer_enable_1,16,onlyOne);
  latch &= compareV1("tkr_out_mask",one.tkr_out_mask,other.tkr_out_mask,16,onlyOne);
 
  latch &= compareV1("conditions",one.conditions,other.conditions,32,onlyOne);
  latch &= compareV1("engine",one.engine,other.engine,16,onlyOne);
  latch &= compareV1("towers",one.towers,other.towers,4,onlyOne);
  latch &= compareV1("tiles",one.tiles,other.tiles,12,onlyOne);
  latch &= compareV1("r",one.r,other.r,54,onlyOne);

  // arrays
  latch &= compareV2("config_reg",(UShort_t*)one.config_reg,(UShort_t*)other.config_reg,12,18,onlyOne);
  latch &= compareV2("tci_dac",(UShort_t*)one.tci_dac,(UShort_t*)other.tci_dac,12,18,onlyOne);

  if ( fullCompare ) {
    std::cout << "Do acd stuff" << std::endl; // sad to make these histos here not in ConfigReport
    TH1* h_pha = new TH1I("delta_pha_threshold","delta_pha_threshold",8192,-4096.,4096.);//binning naively doubled and rounded
    TH1* hv_pha = new TH1S("delta_pha_threshold_v","delta_pha_threshold_v",216,-0.5,215.5);
    l.Add(h_pha); l.Add(hv_pha); 

    TH1* h_bias = new TH1I("delta_bias_dac","delta_bias_dac",16,-8.,8.);
    TH1* hv_bias = new TH1S("delta_bias_dac_v","delta_bias_dac_v",216,-0.5,215.5);
    l.Add(h_bias); l.Add(hv_bias); 

    TH1* h_hld = new TH1I("delta_hld_dac","delta_hld_dac",128,-64.,64.);
    TH1* hv_hld = new TH1S("delta_hld_dac_v","delta_hld_dac_v",512,-216.,216.);
    l.Add(h_hld); l.Add(hv_hld); 

    TH1* h_veto = new TH1I("delta_veto_dac","delta_veto_dac",1024,-63.5,63.5);
    TH1* hv_veto = new TH1F("delta_veto_dac_v","delta_veto_dac_v",512,-216.,216.);
    l.Add(h_veto); l.Add(hv_veto); 

    latch &= compareV2("pha_threshold",(UShort_t*)one.pha_threshold,(UShort_t*)other.pha_threshold,12,18,h_pha, hv_pha, onlyOne);
    latch &= compareV2("bias_dac",(UShort_t*)one.bias_dac,(UShort_t*)other.bias_dac,12,18,h_bias, hv_bias, onlyOne);
    latch &= compareV2("hld_dac",(UShort_t*)one.hld_dac,(UShort_t*)other.hld_dac,12,18,h_hld, hv_hld, onlyOne);
    latch &= compareV2("veto_dac",(UShort_t*)one.veto_dac,(UShort_t*)other.veto_dac,12, 18, h_veto, hv_veto,onlyOne);   
    latch &= compareV2("veto_vernier",(UShort_t*)one.veto_vernier,(UShort_t*)other.veto_vernier,12,18,onlyOne);
  }

  latch &= compareV2("low",(UInt_t*)one.low,(UInt_t*)other.low,16,36,onlyOne);   
  latch &= compareV2("high",(UInt_t*)one.high,(UInt_t*)other.high,16,36,onlyOne); 
  
  latch &= compareV2("tcc_trg_align",(UInt_t*)one.tcc_trg_align,(UInt_t*)other.tcc_trg_align,16,8,onlyOne);   
  latch &= compareV2("tcc_configuration",(UInt_t*)one.tcc_configuration,(UInt_t*)other.tcc_configuration,16,8,onlyOne);   
  latch &= compareV2("input_mask",(UInt_t*)one.input_mask,(UInt_t*)other.input_mask,16,8,onlyOne);  
 
  latch &= compareV2("ccc_trg_alignment",(UInt_t*)one.ccc_trg_alignment,(UInt_t*)other.ccc_trg_alignment,16,4,onlyOne);   
  latch &= compareV2("ccc_configuration",(UInt_t*)one.ccc_configuration,(UInt_t*)other.ccc_configuration,16,4,onlyOne);   
  latch &= compareV2("layer_mask_0",(UInt_t*)one.layer_mask_0,(UInt_t*)other.layer_mask_0,16,4,onlyOne);   
  latch &= compareV2("layer_mask_1",(UInt_t*)one.layer_mask_1,(UInt_t*)other.layer_mask_1,16,4,onlyOne);   

  latch &= compareV3("config",(UInt_t*)one.config,(UInt_t*)other.config,16,4,4,onlyOne);   
  latch &= compareV3("crc_dac",(UInt_t*)one.crc_dac,(UInt_t*)other.crc_dac,16,4,4,onlyOne);   
  latch &= compareV3("delay_1",(UInt_t*)one.delay_1,(UInt_t*)other.delay_1,16,4,4,onlyOne);   
  latch &= compareV3("delay_2",(UInt_t*)one.delay_2,(UInt_t*)other.delay_2,16,4,4,onlyOne);   
  latch &= compareV3("delay_3",(UInt_t*)one.delay_3,(UInt_t*)other.delay_3,16,4,4,onlyOne);   
  
  latch &= compareV3("injection",(UInt_t*)one.injection,(UInt_t*)other.injection,16,36,24,onlyOne);
  latch &= compareV3("calib_mask",(ULong64_t*)one.calib_mask,(ULong64_t*)other.calib_mask,16,36,24,onlyOne);
  
  if ( fullCompare ) {  
    latch &= compareV3("trig_enable",(ULong64_t*)one.trig_enable,(ULong64_t*)other.trig_enable,16,36,24,onlyOne);
    latch &= compareV3("data_mask",(ULong64_t*)one.data_mask,(ULong64_t*)other.data_mask,16,36,24,onlyOne);
    latch &= compareV3("threshold",(UInt_t*)one.threshold,(UInt_t*)other.threshold,16,36,24,onlyOne);
  }

  latch &= compareV3("trc_csr",(ULong64_t*)one.trc_csr,(ULong64_t*)other.trc_csr,16,8,9,onlyOne);
  
  latch &= compareV4("config_0",(UInt_t*)one.config_0,(UInt_t*)other.config_0,16,4,4,12,onlyOne);
  latch &= compareV4("config_1",(UInt_t*)one.config_1,(UInt_t*)other.config_1,16,4,4,12,onlyOne);
  latch &= compareV4("ref_dac",(UInt_t*)one.ref_dac,(UInt_t*)other.ref_dac,16,4,4,12,onlyOne);
  
  if ( fullCompare ) {  
    latch &= compareV4("log_acpt",(UInt_t*)one.log_acpt,(UInt_t*)other.log_acpt,16,4,4,12,onlyOne);
    latch &= compareV4("fle_dac",(UInt_t*)one.fle_dac,(UInt_t*)other.fle_dac,16,4,4,12,onlyOne);
    latch &= compareV4("fhe_dac",(UInt_t*)one.fhe_dac,(UInt_t*)other.fhe_dac,16,4,4,12,onlyOne);
    latch &= compareV4("rng_uld_dac",(UInt_t*)one.rng_uld_dac,(UInt_t*)other.rng_uld_dac,16,4,4,12,onlyOne);
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
    printf("%s 0x%0x 0x%0x\n",name,(UShort_t)val1,(UShort_t)val2);
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

Bool_t ConfigCompare::compareV2(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, TH1* h, TH1* hv, Bool_t onlyOne) {
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      int idx = i*size2 + j;
      if ( val1[idx] != val2[idx] ) {
	printf("%s[%i,%i] 0x%0x 0x%0x\n",name,i,j,val1[idx],val2[idx]);
	//	if (onlyOne) return kFALSE;
	h->Fill(val1[idx]-val2[idx]);
	hv->SetBinContent(idx+1,val1[idx]-val2[idx]);
      }
    }
  }
  return kTRUE;
}

Bool_t ConfigCompare::compareV2(const char* name, UShort_t* val1, UShort_t* val2, int size1, int size2, TH1* h, TH1* hv, Bool_t onlyOne){
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      int idx = i*size2 + j;
      if ( val1[idx] != val2[idx] ) {
	printf("%s[%i,%i] 0x%0x 0x%0x\n",name,i,j,val1[idx],val2[idx]);
       	//if (onlyOne) return kFALSE;
	h->Fill(val1[idx]-val2[idx]);
	hv->SetBinContent(idx+1,val1[idx]-val2[idx]);
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
  


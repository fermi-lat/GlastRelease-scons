#define ConfigReport_cxx
#include "configData/ConfigReport.h"
#include "configData/ConfigBase.h"

#include <TH2.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH1S.h>
#include <TList.h>
#include <stdio.h>
#include <iostream>


void     
ConfigReport::report(ConfigBase& c, TList& l, Option_t* options, std::ostream& os) {

  reportSingleton("aem_configuration",c.aem_configuration,os);
  reportSingleton("trgseq",c.trgseq,os);
  reportSingleton("window_width",c.window_width,os);
  reportSingleton("window_open_mask",c.window_open_mask,os);
  reportSingleton("periodic_prescale",c.periodic_prescale,os);
  reportSingleton("periodic_limit",c.periodic_limit,os);
  reportSingleton("periodic_free",c.periodic_free,os);
  reportSingleton("periodic_onepps",c.periodic_onepps,os);
  reportSingleton("acd_cno",c.acd_cno,os);
  reportSingleton("tower_busy",c.tower_busy,os);
  reportSingleton("external",c.external,os);
  reportSingleton("configuration",c.configuration,os);

  reportV1("hitmap_delay",c.hitmap_delay_bcast,c.hitmap_delay,12,os);  
  reportV1("hitmap_width",c.hitmap_width_bcast,c.hitmap_width,12,os);
  reportV1("hitmap_deadtime",c.hitmap_deadtime_bcast,c.hitmap_deadtime,12,os);
  reportV1("veto_delay",c.veto_delay_bcast,c.veto_delay,12,os);
  reportV1("veto_width",c.veto_width_bcast,c.veto_width,12,os);
  reportV1("adc_tacq",c.adc_tacq_bcast,c.adc_tacq,12,os);
  reportV1("data_masks",c.data_masks_bcast,c.data_masks,16,os);
  reportV1("tkr_trgseq",c.tkr_trgseq_bcast,c.tkr_trgseq,16,os);
  reportV1("cal_trgseq",c.cal_trgseq_bcast,c.cal_trgseq,16,os); 

  reportV1("cal_in_mask",c.cal_in_mask_bcast,c.cal_in_mask,16,os); 
  reportV1("tkr_layer_enable_0",c.tkr_layer_enable_0_bcast,c.tkr_layer_enable_0,16,os); 
  reportV1("tkr_layer_enable_1",c.tkr_layer_enable_1_bcast,c.tkr_layer_enable_1,16,os); 
  reportV1("tkr_out_mask",c.tkr_out_mask_bcast,c.tkr_out_mask,16,os); 

  //reportV1("conditions",c.conditions_bcast,c.conditions,32,os);
  //reportV1("engine",c.engine_bcast,c.engine,16,os);
  //reportV1("towers",c.towers_bcast,c.towers,4,os);
  //reportV1("tiles",c.tiles_bcast,c.tiles,12,os);
  //reportV1("r",c.r_bcast,c.r,54,os);

  reportV2("tci_dac",c.tci_dac_bcast,(UShort_t*)c.tci_dac,12,18,os);
  reportV2("config_reg",c.config_reg_bcast,(UShort_t*)c.config_reg,12,18,os);

  TH1* h_pha = new TH1I("pha_threshold","pha_threshold",4096,-0.5,4095.5);
  TH1* hv_pha = new TH1S("pha_threshold_v","pha_threshold_v",216,-0.5,215.5);
  l.Add(h_pha); l.Add(hv_pha); 
  
  TH1* h_bias = new TH1I("bias_dac","bias_dac",8,-0.5,7.5);
  TH1* hv_bias = new TH1S("bias_dac_v","bias_dac_v",216,-0.5,215.5);
  l.Add(h_bias); l.Add(hv_bias); 

  TH1* h_hld = new TH1I("hld_dac","hld_dac",64,-0.5,63.5);
  TH1* hv_hld = new TH1S("hld_dac_v","hld_dac_v",216,-0.5,215.5);
  l.Add(h_hld); l.Add(hv_hld); 

  TH1* h_veto = new TH1I("veto_dac","veto_dac",1024,-0.5,63.5);
  TH1* hv_veto = new TH1F("veto_dac_v","veto_dac_v",216,-0.5,215.5);
  l.Add(h_veto); l.Add(hv_veto); 

  reportV2("pha_threshold",c.pha_threshold_bcast,(UShort_t*)c.pha_threshold,12,18,os,h_pha,hv_pha);
  reportV2("bias_dac",c.bias_dac_bcast,(UShort_t*)c.bias_dac,12,18,os,h_bias,hv_bias);
  reportV2("hld_dac",c.hld_dac_bcast,(UShort_t*)c.hld_dac,12,18,os,h_hld,hv_hld);
  reportVeto("veto_dac",(UShort_t*)c.veto_dac,(UShort_t*)c.veto_vernier,os,h_veto,hv_veto);   

  reportV2("low",c.low_bcast,(UInt_t*)c.low,16,36,os);   
  reportV2("high",c.high_bcast,(UInt_t*)c.high,16,36,os); 
  
  reportV2("tcc_trg_align",c.tcc_trg_align_bcast,(UInt_t*)c.tcc_trg_align,16,8,os);   
  reportV2("tcc_configuration",c.tcc_configuration_bcast,(UInt_t*)c.tcc_configuration,16,8,os);   
  reportV2("input_mask",c.input_mask_bcast,(UInt_t*)c.input_mask,16,8,os);  
 
  reportV2("ccc_trg_alignment",c.ccc_trg_alignment_bcast,(UInt_t*)c.ccc_trg_alignment,16,4,os);   
  reportV2("ccc_configuration",c.ccc_configuration_bcast,(UInt_t*)c.ccc_configuration,16,4,os);   
  reportV2("layer_mask_0",c.layer_mask_0_bcast,(UInt_t*)c.layer_mask_0,16,4,os);   
  reportV2("layer_mask_1",c.layer_mask_1_bcast,(UInt_t*)c.layer_mask_1,16,4,os);   

  reportV2("crc_dac",c.crc_dac_bcast,(UInt_t*)c.crc_dac,16,4,os);   
  reportV2("config",c.config_bcast,(UInt_t*)c.config,16,4,os);   
  reportV2("delay_1",c.delay_1_bcast,(UInt_t*)c.delay_1,16,4,os);   
  reportV2("delay_2",c.delay_2_bcast,(UInt_t*)c.delay_2,16,4,os);   
  reportV2("delay_3",c.delay_3_bcast,(UInt_t*)c.delay_3,16,4,os);   

  TH1* h_trig_enable = new TH1I("trig_enable","trig_enable",64,-0.5,64.5);
  TH1* hv_trig_enable = new TH1S("trig_enable_v","trig_enable_v",13824,-0.5,13823.5);  
  TH1* h_data_mask = new TH1I("data_mask","data_mask",64,-0.5,64.5);
  TH1* hv_data_mask = new TH1S("data_mask_v","data_mask_v",13824,-0.5,13823.5);  
  reportMask("trig_enable",c.trig_enable_bcast,(ULong64_t*)c.trig_enable,16,36,24,os,h_trig_enable,hv_trig_enable);
  reportMask("data_mask",c.data_mask_bcast,(ULong64_t*)c.data_mask,16,36,24,os,h_data_mask,hv_data_mask);
  reportMask("calib_mask",c.calib_mask_bcast,(ULong64_t*)c.calib_mask,16,36,24,os);
  l.Add(h_trig_enable); l.Add(hv_trig_enable); 
  l.Add(h_data_mask); l.Add(hv_data_mask);   

  reportV3("injection",c.injection_bcast,(UInt_t*)c.injection,16,36,24,os);
  TH1* h_thresh = new TH1I("threshold","threshold",128,-0.5,127.5);
  TH1* hv_thresh = new TH1S("threshold_v","threshold_v",13824,-0.5,13823.5);
  l.Add(h_thresh); l.Add(hv_thresh); 

  reportV3("threshold",c.threshold_bcast,(UInt_t*)c.threshold,16,36,24,os,h_thresh,hv_thresh);
  reportV3("trc_csr",c.trc_csr_bcast,(ULong64_t*)c.trc_csr,16,8,9,os);
  
  reportV4("config_0",c.config_0_bcast,(UInt_t*)c.config_0,16,4,4,12,os);
  reportV4("config_1",c.config_1_bcast,(UInt_t*)c.config_1,16,4,4,12,os);
  reportV4("ref_dac",c.ref_dac_bcast,(UInt_t*)c.ref_dac,16,4,4,12,os);

  TH1* h_log_acpt = new TH1I("log_acpt","log_acpt",128,-0.5,127.5);
  TH1* hv_log_acpt = new TH1S("log_acpt_v","log_acpt_v",3072,-0.5,3071.5);
  l.Add(h_log_acpt); l.Add(hv_log_acpt); 
  TH1* h_fle = new TH1I("fle","fle",128,-0.5,127.5);
  TH1* hv_fle = new TH1S("fle_v","fle_v",3072,-0.5,3071.5);
  l.Add(h_fle); l.Add(hv_fle); 
  TH1* h_fhe = new TH1I("fhe","fhe",128,-0.5,127.5);
  TH1* hv_fhe = new TH1S("fhe_v","fhe_v",3072,-0.5,3071.5);
  l.Add(h_fhe); l.Add(hv_fhe); 
  TH1* h_rng_uld = new TH1I("rng_uld","rng_uld",128,-0.5,127.5);
  TH1* hv_rng_uld = new TH1S("rng_uld_v","rng_uld_v",3072,-0.5,3071.5); 
  l.Add(h_rng_uld); l.Add(hv_rng_uld); 
  
  reportV4("log_acpt",c.log_acpt_bcast,(UInt_t*)c.log_acpt,16,4,4,12,os,h_log_acpt,hv_log_acpt);
  reportV4("fle_dac",c.fle_dac_bcast,(UInt_t*)c.fle_dac,16,4,4,12,os,h_fle,hv_fle);
  reportV4("fhe_dac",c.fhe_dac_bcast,(UInt_t*)c.fhe_dac,16,4,4,12,os,h_fhe,hv_fhe);
  reportV4("rng_uld_dac",c.rng_uld_dac_bcast,(UInt_t*)c.rng_uld_dac,16,4,4,12,os,h_rng_uld,hv_rng_uld);

}


void ConfigReport::reportSingleton(const char* name, ULong64_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << val << std::endl;
}

void ConfigReport::reportSingleton(const char* name, UInt_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << val << std::endl;
}

void ConfigReport::reportSingleton(const char* name, Char_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << (UShort_t)val << std::endl;
}

void ConfigReport::reportV1(const char* name, UInt_t bcast, UInt_t* val, int size, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  for (int i(0); i < size; i++ ) {
    if ( h == 0 ) {
      if ( val[i] != bcast ) {
	printf("%s[%i] 0x%0x\n",name,i,val[i]);
      }
    } else {
      h->Fill(val[i]);
      hv->SetBinContent(i+1,val[i]);
    }
  }
}

void ConfigReport::reportV1(const char* name, UShort_t bcast, UShort_t* val, int size, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  for (int i(0); i < size; i++ ) {
    if ( h == 0 ) {
      if ( val[i] != bcast ) {
	printf("%s[%i] 0x%0x\n",name,i,val[i]);
      }
    } else {
      h->Fill(val[i]);
      hv->SetBinContent(i+1,val[i]);
    }
  }
}

void ConfigReport::reportV2(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, 
			     std::ostream& os,  TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      if ( h == 0 ) {
	if ( val[idx] != bcast ) {
	  printf("%s[%i,%i] 0x%0x\n",name,i,j,val[idx]);
	}
      } else {
	h->Fill(val[idx]);
	hv->SetBinContent(idx+1,val[idx]);
      }
      idx++;
    }
  }
}

void ConfigReport::reportV2(const char* name, UShort_t bcast, UShort_t* val, int size1, int size2, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      if ( h == 0 ) {
	if ( val[idx] != bcast ) {
	  printf("%s[%i,%i] 0x%0x\n",name,i,j,val[idx]);
	}
      } else {
	h->Fill(val[idx]);
	hv->SetBinContent(idx+1,val[idx]);
      }
      idx++;
    }
  }
}

void ConfigReport::reportVeto(const char* name, UShort_t* veto, UShort_t* vernier,
			       std::ostream& os, TH1* h, TH1* hv) {
  Int_t idx(0);
  for (int i(0); i < 12; i++ ) {
    for (int j(0); j < 18; j++ ) {
      Float_t vx = (Float_t)(veto[idx]) + (((Float_t)(vernier[idx]))/32.);
      h->Fill(vx);
      hv->SetBinContent(idx+1,vx);
      idx++;
    }
  }
}

void ConfigReport::reportV3(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {      
	if ( h == 0 ) {
	  if ( val[idx] != bcast ) {
	    printf("%s[%i,%i,%i] 0x%0x\n",name,i,j,k,val[idx]);
	  }
	} else {
	  h->Fill(val[idx]);
	  hv->SetBinContent(idx+1,val[idx]);
	}
	idx++;
      }
    }
  }
}

void ConfigReport::reportV3(const char* name, ULong64_t bcast, ULong64_t* val, int size1, int size2, int size3, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {      
	if ( h == 0 ) {
	  if ( val[idx] != bcast ) {
	    printf("%s[%i,%i,%i] 0x%0x\n",name,i,j,k,val[idx]);
	  }
	} else {
	  h->Fill(val[idx]);
	  hv->SetBinContent(idx+1,val[idx]);
	}
	idx++;
      }
    }
  }
}

void ConfigReport::reportMask(const char* name, ULong64_t bcast, ULong64_t* val, int size1, int size2, int size3, 
			       std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {      
	if ( h == 0 ) {
	  if ( val[idx] != bcast ) {
	    os << name << '[' << i << ',' << j << ',' << k << ']' << std::hex << val[idx] << std::endl;
	  }
	} else {
	  Int_t nDiff(0);
	  Long64_t theMask(0x1);
	  for ( int bit(0); bit < 64; bit++ ) {
	    if ( ( bcast & theMask ) != ( val[idx] & theMask ) ) {
	      nDiff++;
	    }
	    theMask = theMask << 1;
	  }
	  h->Fill(nDiff);
	  hv->SetBinContent(idx+1,nDiff);
	}
	idx++;
      }
    }
  }
}

void ConfigReport::reportV4(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, int size4, 
			     std::ostream& os, TH1* h, TH1* hv) {
  os << "\t" << name << " BCAST: 0x" << std::hex << bcast << std::endl;
  Int_t idx(0);
  for (int i(0); i < size1; i++ ) {
    for (int j(0); j < size2; j++ ) {
      for (int k(0); k < size3; k++ ) {      
	for (int l(0); l < size4; l++ ) {      
	  if ( h == 0 ) {
	    if ( val[idx] != bcast ) {
	      printf("%s[%i,%i,%i,%i] 0x%0x\n",name,i,j,k,l,val[idx]);
	    }
	  } else {
	    h->Fill(val[idx]);
	    hv->SetBinContent(idx+1,val[idx]);
	  }
	  idx++;
	}
      }
    }
  }
}



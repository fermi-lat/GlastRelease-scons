#define ConfigCompare_cxx
#include "configData/ConfigCompare.h"
#include <TH2.h>
#include <TH1.h>
#include <TH1I.h>
#include <TH1S.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TList.h>
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

void     
ConfigCompare::report(TList& l, Option_t* options, std::ostream& os) {

  LoadTree(0); GetEntry(0);

  reportSingleton("aem_configuration",aem_configuration,os);
  reportSingleton("trgseq",trgseq,os);
  reportSingleton("window_width",window_width,os);
  reportSingleton("window_open_mask",window_open_mask,os);
  reportSingleton("periodic_prescale",periodic_prescale,os);
  reportSingleton("periodic_limit",periodic_limit,os);
  reportSingleton("periodic_free",periodic_free,os);
  reportSingleton("periodic_onepps",periodic_onepps,os);
  reportSingleton("acd_cno",acd_cno,os);
  reportSingleton("tower_busy",tower_busy,os);
  reportSingleton("external",external,os);
  reportSingleton("configuration",configuration,os);

  reportV1("hitmap_delay",hitmap_delay_bcast,hitmap_delay,12,os);  
  reportV1("hitmap_width",hitmap_width_bcast,hitmap_width,12,os);
  reportV1("hitmap_deadtime",hitmap_deadtime_bcast,hitmap_deadtime,12,os);
  reportV1("veto_delay",veto_delay_bcast,veto_delay,12,os);
  reportV1("veto_width",veto_width_bcast,veto_width,12,os);
  reportV1("adc_tacq",adc_tacq_bcast,adc_tacq,12,os);
  reportV1("data_masks",data_masks_bcast,data_masks,16,os);
  reportV1("tkr_trgseq",tkr_trgseq_bcast,tkr_trgseq,16,os);
  reportV1("cal_trgseq",cal_trgseq_bcast,cal_trgseq,16,os); 

  reportV1("cal_in_mask",cal_in_mask_bcast,cal_in_mask,16,os); 
  reportV1("tkr_layer_enable_0",tkr_layer_enable_0_bcast,tkr_layer_enable_0,16,os); 
  reportV1("tkr_layer_enable_1",tkr_layer_enable_1_bcast,tkr_layer_enable_1,16,os); 
  reportV1("tkr_out_mask",tkr_out_mask_bcast,tkr_out_mask,16,os); 

  //reportV1("conditions",conditions_bcast,conditions,32,os);
  //reportV1("engine",engine_bcast,engine,16,os);
  //reportV1("towers",towers_bcast,towers,4,os);
  //reportV1("tiles",tiles_bcast,tiles,12,os);
  //reportV1("r",r_bcast,r,54,os);

  reportV2("tci_dac",tci_dac_bcast,(UShort_t*)tci_dac,12,18,os);
  reportV2("config_reg",config_reg_bcast,(UShort_t*)config_reg,12,18,os);

  TH1* h_pha = new TH1I("pha_threshold","pha_threshold",4096,-0.5,4095.5);
  TH1* hv_pha = new TH1S("pha_threshold_v","pha_threshold_v",216,-0.5,215.5);
  l.Add(h_pha); l.Add(hv_pha); 
  
  TH1* h_bias = new TH1I("bias_dac","bias_dac",8,-0.5,7.5);
  TH1* hv_bias = new TH1S("bias_dac_v","bias_dac_v",216,-0.5,215.5);
  l.Add(h_bias); l.Add(hv_bias); 

  TH1* h_hld = new TH1I("hld_dac","hld_dac",64,-0.5,63.5);
  TH1* hv_hld = new TH1S("hld_dac_v","hld_dac_v",216,-0.5,215.5);
  l.Add(h_bias); l.Add(hv_bias); 

  TH1* h_veto = new TH1I("veto_dac","veto_dac",1024,-0.5,63.5);
  TH1* hv_veto = new TH1F("veto_dac_v","veto_dac_v",216,-0.5,215.5);
  l.Add(h_veto); l.Add(hv_veto); 

  reportV2("pha_threshold",pha_threshold_bcast,(UShort_t*)pha_threshold,12,18,os,h_pha,hv_pha);
  reportV2("bias_dac",bias_dac_bcast,(UShort_t*)bias_dac,12,18,os,h_bias,hv_bias);
  reportV2("hld_dac",hld_dac_bcast,(UShort_t*)hld_dac,12,18,os,h_hld,hv_hld);
  reportVeto("veto_dac",(UShort_t*)veto_dac,(UShort_t*)veto_vernier,os,h_veto,hv_veto);   

  reportV2("low",low_bcast,(UInt_t*)low,16,36,os);   
  reportV2("high",high_bcast,(UInt_t*)high,16,36,os); 
  
  reportV2("tcc_trg_align",tcc_trg_align_bcast,(UInt_t*)tcc_trg_align,16,8,os);   
  reportV2("tcc_configuration",tcc_configuration_bcast,(UInt_t*)tcc_configuration,16,8,os);   
  reportV2("input_mask",input_mask_bcast,(UInt_t*)input_mask,16,8,os);  
 
  reportV2("ccc_trg_alignment",ccc_trg_alignment_bcast,(UInt_t*)ccc_trg_alignment,16,4,os);   
  reportV2("ccc_configuration",ccc_configuration_bcast,(UInt_t*)ccc_configuration,16,4,os);   
  reportV2("layer_mask_0",layer_mask_0_bcast,(UInt_t*)layer_mask_0,16,4,os);   
  reportV2("layer_mask_1",layer_mask_1_bcast,(UInt_t*)layer_mask_1,16,4,os);   

  reportV2("crc_dac",crc_dac_bcast,(UInt_t*)crc_dac,16,4,os);   
  reportV2("config",config_bcast,(UInt_t*)config,16,4,os);   
  reportV2("delay_1",delay_1_bcast,(UInt_t*)delay_1,16,4,os);   
  reportV2("delay_2",delay_2_bcast,(UInt_t*)delay_2,16,4,os);   
  reportV2("delay_3",delay_3_bcast,(UInt_t*)delay_3,16,4,os);   

  TH1* h_trig_enable = new TH1I("trig_enable","trig_enable",64,-0.5,64.5);
  TH1* hv_trig_enable = new TH1S("trig_enable_v","trig_enable_v",13824,-0.5,13823.5);  
  TH1* h_data_mask = new TH1I("data_mask","data_mask",64,-0.5,64.5);
  TH1* hv_data_mask = new TH1S("data_mask_v","data_mask_v",13824,-0.5,13823.5);  
  reportMask("trig_enable",trig_enable_bcast,(ULong64_t*)trig_enable,16,36,24,os,h_trig_enable,hv_trig_enable);
  reportMask("data_mask",data_mask_bcast,(ULong64_t*)data_mask,16,36,24,os,h_data_mask,hv_data_mask);
  reportMask("calib_mask",calib_mask_bcast,(ULong64_t*)calib_mask,16,36,24,os);
  l.Add(h_trig_enable); l.Add(hv_trig_enable); 
  l.Add(h_data_mask); l.Add(hv_data_mask);   

  reportV3("injection",injection_bcast,(UInt_t*)injection,16,36,24,os);
  TH1* h_thresh = new TH1I("threshold","threshold",128,-0.5,127.5);
  TH1* hv_thresh = new TH1S("threshold_v","threshold_v",13824,-0.5,13823.5);
  l.Add(h_thresh); l.Add(hv_thresh); 

  reportV3("threshold",threshold_bcast,(UInt_t*)threshold,16,36,24,os,h_thresh,hv_thresh);
  reportV3("trc_csr",trc_csr_bcast,(UInt_t*)trc_csr,16,8,9,os);
  
  reportV4("config_0",config_0_bcast,(UInt_t*)config_0,16,4,4,12,os);
  reportV4("config_1",config_1_bcast,(UInt_t*)config_1,16,4,4,12,os);
  reportV4("ref_dac",ref_dac_bcast,(UInt_t*)ref_dac,16,4,4,12,os);

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
  
  reportV4("log_acpt",log_acpt_bcast,(UInt_t*)log_acpt,16,4,4,12,os,h_log_acpt,hv_log_acpt);
  reportV4("fle_dac",fle_dac_bcast,(UInt_t*)fle_dac,16,4,4,12,os,h_fle,hv_fle);
  reportV4("fhe_dac",fhe_dac_bcast,(UInt_t*)fhe_dac,16,4,4,12,os,h_fhe,hv_fhe);
  reportV4("rng_uld_dac",rng_uld_dac_bcast,(UInt_t*)rng_uld_dac,16,4,4,12,os,h_rng_uld,hv_rng_uld);

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
  latch &= compareV1("hold_delay",hold_delay,other.hold_delay,12);
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
  

void ConfigCompare::reportSingleton(const char* name, ULong64_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << val << std::endl;
}

void ConfigCompare::reportSingleton(const char* name, UInt_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << val << std::endl;
}

void ConfigCompare::reportSingleton(const char* name, Char_t val, std::ostream& os) {
  os << "\t" << name << " : 0x" << std::hex << (UShort_t)val << std::endl;
}

void ConfigCompare::reportV1(const char* name, UInt_t bcast, UInt_t* val, int size, 
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

void ConfigCompare::reportV1(const char* name, UShort_t bcast, UShort_t* val, int size, 
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

void ConfigCompare::reportV2(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, 
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

void ConfigCompare::reportV2(const char* name, UShort_t bcast, UShort_t* val, int size1, int size2, 
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

void ConfigCompare::reportVeto(const char* name, UShort_t* veto, UShort_t* vernier,
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

void ConfigCompare::reportV3(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, 
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

void ConfigCompare::reportMask(const char* name, ULong64_t bcast, ULong64_t* val, int size1, int size2, int size3, 
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

void ConfigCompare::reportV4(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, int size4, 
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
  


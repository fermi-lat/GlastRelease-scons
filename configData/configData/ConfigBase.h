//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 25 15:10:52 2007 by ROOT version 5.10/00
// from TTree ConfigBase/ConfigBaseuration Tree
// found on file: e2e_LAT-20xGammaThruMIPCNONoPer_0.01hr.root
//////////////////////////////////////////////////////////

#ifndef ConfigBase_h
#define ConfigBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>

class ConfigBase {
public :

   ConfigBase(TTree *tree=0);
   virtual ~ConfigBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   UInt_t          aem_configuration;
   UInt_t          trgseq;
   UShort_t        hold_delay_bcast;
   UShort_t        hold_delay[12];
   UShort_t        hitmap_delay_bcast;
   UShort_t        hitmap_delay[12];
   UShort_t        hitmap_width_bcast;
   UShort_t        hitmap_width[12];
   UShort_t        hitmap_deadtime_bcast;
   UShort_t        hitmap_deadtime[12];
   UShort_t        veto_delay_bcast;
   UShort_t        veto_delay[12];
   UShort_t        veto_width_bcast;
   UShort_t        veto_width[12];
   UShort_t        adc_tacq_bcast;
   UShort_t        adc_tacq[12];
   UShort_t        pha_threshold_bcast;
   UShort_t        pha_threshold[12][18];
   UShort_t        config_reg_bcast;
   UShort_t        config_reg[12][18];
   UShort_t        tci_dac_bcast;
   UShort_t        tci_dac[12][18];
   UShort_t        bias_dac_bcast;
   UShort_t        bias_dac[12][18];
   UShort_t        hld_dac_bcast;
   UShort_t        hld_dac[12][18];
   UShort_t        veto_dac_bcast;
   UShort_t        veto_dac[12][18];
   UShort_t        veto_vernier_bcast;
   UShort_t        veto_vernier[12][18];
   UInt_t          data_masks_bcast;
   UInt_t          data_masks[16];
   UInt_t          tkr_trgseq_bcast;
   UInt_t          tkr_trgseq[16];
   UInt_t          cal_trgseq_bcast;
   UInt_t          cal_trgseq[16];   
   UInt_t          cal_in_mask_bcast;
   UInt_t          cal_in_mask[16];
   UInt_t          tkr_layer_enable_0_bcast;
   UInt_t          tkr_layer_enable_0[16];
   UInt_t          tkr_layer_enable_1_bcast;
   UInt_t          tkr_layer_enable_1[16];
   UInt_t          tkr_out_mask_bcast;
   UInt_t          tkr_out_mask[16];
   UInt_t          low_bcast;
   UInt_t          low[16][36];
   UInt_t          high_bcast;
   UInt_t          high[16][36];
   ULong64_t       trig_enable_bcast;
   ULong64_t       trig_enable[16][36][24];
   ULong64_t       data_mask_bcast;
   ULong64_t       data_mask[16][36][24];
   ULong64_t       calib_mask_bcast;
   ULong64_t       calib_mask[16][36][24];
   UInt_t          threshold_bcast;
   UInt_t          threshold[16][36][24];
   UInt_t          injection_bcast;
   UInt_t          injection[16][36][24];
   UInt_t          tcc_trg_align_bcast;
   UInt_t          tcc_trg_align[16][8];
   UInt_t          tcc_configuration_bcast;
   UInt_t          tcc_configuration[16][8];
   UInt_t          input_mask_bcast;
   UInt_t          input_mask[16][8];
   ULong64_t       trc_csr_bcast;
   ULong64_t       trc_csr[16][8][9];
   UInt_t          ccc_trg_alignment_bcast;
   UInt_t          ccc_trg_alignment[16][4];
   UInt_t          ccc_configuration_bcast;
   UInt_t          ccc_configuration[16][4];
   UInt_t          layer_mask_0_bcast;
   UInt_t          layer_mask_0[16][4];
   UInt_t          layer_mask_1_bcast;
   UInt_t          layer_mask_1[16][4];
   UInt_t          config_bcast;
   UInt_t          config[16][4][4];
   UInt_t          crc_dac_bcast;
   UInt_t          crc_dac[16][4][4];
   UInt_t          delay_1_bcast;
   UInt_t          delay_1[16][4][4];
   UInt_t          delay_2_bcast;
   UInt_t          delay_2[16][4][4];
   UInt_t          delay_3_bcast;
   UInt_t          delay_3[16][4][4];
   UInt_t          config_0_bcast;
   UInt_t          config_0[16][4][4][12];
   UInt_t          config_1_bcast;
   UInt_t          config_1[16][4][4][12];
   UInt_t          ref_dac_bcast;
   UInt_t          ref_dac[16][4][4][12];
   UInt_t          log_acpt_bcast;
   UInt_t          log_acpt[16][4][4][12];
   UInt_t          fle_dac_bcast;
   UInt_t          fle_dac[16][4][4][12];
   UInt_t          fhe_dac_bcast;
   UInt_t          fhe_dac[16][4][4][12];
   UInt_t          rng_uld_dac_bcast;
   UInt_t          rng_uld_dac[16][4][4][12];
   UInt_t          window_width;
   Char_t          window_open_mask;
   UInt_t          periodic_prescale;
   UInt_t          periodic_limit;
   Char_t          periodic_free;
   Char_t          periodic_onepps;
   UInt_t          conditions[32];
   UInt_t          engine[16];
   UInt_t          towers[4];
   UInt_t          tiles[12];
   UInt_t          acd_cno;
   UInt_t          tower_busy;
   Char_t          external;
   UInt_t          configuration;
   UInt_t          r[54];

protected:
  
private:

   // List of branches
   TBranch        *b_aem_configuration;   //!
   TBranch        *b_trgseq;   //!

   TBranch        *b_hold_delay_bcast;   //!
   TBranch        *b_hitmap_delay_bcast;   //!
   TBranch        *b_hitmap_width_bcast;   //!
   TBranch        *b_hitmap_deadtime_bcast;   //!
   TBranch        *b_veto_delay_bcast;   //!
   TBranch        *b_veto_width_bcast;   //!
   TBranch        *b_adc_tacq_bcast;   //!
   TBranch        *b_pha_threshold_bcast;   //!
   TBranch        *b_config_reg_bcast;   //!
   TBranch        *b_tci_dac_bcast;   //!
   TBranch        *b_bias_dac_bcast;   //!
   TBranch        *b_hld_dac_bcast;   //!
   TBranch        *b_veto_dac_bcast;   //!
   TBranch        *b_veto_vernier_bcast;   //!
   TBranch        *b_data_masks_bcast;   //!
   TBranch        *b_tkr_trgseq_bcast;   //!
   TBranch        *b_cal_trgseq_bcast;   //!
   TBranch        *b_cal_in_mask_bcast;   //!
   TBranch        *b_tkr_layer_enable_0_bcast;   //!
   TBranch        *b_tkr_layer_enable_1_bcast;   //!
   TBranch        *b_tkr_out_mask_bcast;   //!
   TBranch        *b_low_bcast;   //!
   TBranch        *b_high_bcast;   //!
   TBranch        *b_trig_enable_bcast;   //!
   TBranch        *b_data_mask_bcast;   //!
   TBranch        *b_calib_mask_bcast;   //!
   TBranch        *b_threshold_bcast;   //!
   TBranch        *b_injection_bcast;   //!
   TBranch        *b_tcc_trg_align_bcast;   //!
   TBranch        *b_tcc_configuration_bcast;   //!
   TBranch        *b_input_mask_bcast;   //!
   TBranch        *b_trc_csr_bcast;   //!
   TBranch        *b_ccc_trg_alignment_bcast;   //!
   TBranch        *b_ccc_configuration_bcast;   //!
   TBranch        *b_layer_mask_0_bcast;   //!
   TBranch        *b_layer_mask_1_bcast;   //!
   TBranch        *b_config_bcast;   //!
   TBranch        *b_crc_dac_bcast;   //!
   TBranch        *b_delay_1_bcast;   //!
   TBranch        *b_delay_2_bcast;   //!
   TBranch        *b_delay_3_bcast;   //!
   TBranch        *b_config_0_bcast;   //!
   TBranch        *b_config_1_bcast;   //!
   TBranch        *b_ref_dac_bcast;   //!
   TBranch        *b_log_acpt_bcast;   //!
   TBranch        *b_fle_dac_bcast;   //!
   TBranch        *b_fhe_dac_bcast;   //!
   TBranch        *b_rng_uld_dac_bcast;   //!


   TBranch        *b_hold_delay;   //!
   TBranch        *b_hitmap_delay;   //!
   TBranch        *b_hitmap_width;   //!
   TBranch        *b_hitmap_deadtime;   //!
   TBranch        *b_veto_delay;   //!
   TBranch        *b_veto_width;   //!
   TBranch        *b_adc_tacq;   //!
   TBranch        *b_pha_threshold;   //!
   TBranch        *b_config_reg;   //!
   TBranch        *b_tci_dac;   //!
   TBranch        *b_bias_dac;   //!
   TBranch        *b_hld_dac;   //!
   TBranch        *b_veto_dac;   //!
   TBranch        *b_veto_vernier;   //!
   TBranch        *b_data_masks;   //!
   TBranch        *b_tkr_trgseq;   //!
   TBranch        *b_cal_trgseq;   //!
   TBranch        *b_cal_in_mask;   //!
   TBranch        *b_tkr_layer_enable_0;   //!
   TBranch        *b_tkr_layer_enable_1;   //!
   TBranch        *b_tkr_out_mask;   //!
   TBranch        *b_low;   //!
   TBranch        *b_high;   //!
   TBranch        *b_trig_enable;   //!
   TBranch        *b_data_mask;   //!
   TBranch        *b_calib_mask;   //!
   TBranch        *b_threshold;   //!
   TBranch        *b_injection;   //!
   TBranch        *b_tcc_trg_align;   //!
   TBranch        *b_tcc_configuration;   //!
   TBranch        *b_input_mask;   //!
   TBranch        *b_trc_csr;   //!
   TBranch        *b_ccc_trg_alignment;   //!
   TBranch        *b_ccc_configuration;   //!
   TBranch        *b_layer_mask_0;   //!
   TBranch        *b_layer_mask_1;   //!
   TBranch        *b_crc_dac;   //!
   TBranch        *b_config;   //!
   TBranch        *b_delay_1;   //!
   TBranch        *b_delay_2;   //!
   TBranch        *b_delay_3;   //!
   TBranch        *b_config_0;   //!
   TBranch        *b_config_1;   //!
   TBranch        *b_ref_dac;   //!
   TBranch        *b_log_acpt;   //!
   TBranch        *b_fle_dac;   //!
   TBranch        *b_fhe_dac;   //!
   TBranch        *b_rng_uld_dac;   //!
   TBranch        *b_window_width;   //!
   TBranch        *b_window_open_mask;   //!
   TBranch        *b_periodic_prescale;   //!
   TBranch        *b_periodic_limit;   //!
   TBranch        *b_periodic_free;   //!
   TBranch        *b_conditions;   //!
   TBranch        *b_periodic_onepps;   //!
   TBranch        *b_engine;   //!
   TBranch        *b_towers;   //!
   TBranch        *b_tiles;   //!
   TBranch        *b_acd_cno;   //!
   TBranch        *b_tower_busy;   //!
   TBranch        *b_external;   //!
   TBranch        *b_configuration;   //!
   TBranch        *b_r;   //!

};

#endif

#ifdef ConfigBase_cxx
ConfigBase::ConfigBase(TTree *tree)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  Init(tree);
}

ConfigBase::~ConfigBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ConfigBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ConfigBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ConfigBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("aem_configuration",&aem_configuration);
   fChain->SetBranchAddress("trgseq",&trgseq);
   fChain->SetBranchAddress("hold_delay",hold_delay);
   fChain->SetBranchAddress("hitmap_delay",hitmap_delay);
   fChain->SetBranchAddress("hitmap_width",hitmap_width);
   fChain->SetBranchAddress("hitmap_deadtime",hitmap_deadtime);
   fChain->SetBranchAddress("veto_delay",veto_delay);
   fChain->SetBranchAddress("veto_width",veto_width);
   fChain->SetBranchAddress("adc_tacq",adc_tacq);
   fChain->SetBranchAddress("pha_threshold",pha_threshold);
   fChain->SetBranchAddress("config_reg",config_reg);
   fChain->SetBranchAddress("tci_dac",tci_dac);
   fChain->SetBranchAddress("bias_dac",bias_dac);
   fChain->SetBranchAddress("hld_dac",hld_dac);
   fChain->SetBranchAddress("veto_dac",veto_dac);
   fChain->SetBranchAddress("veto_vernier",veto_vernier);
   fChain->SetBranchAddress("data_masks",data_masks);
   fChain->SetBranchAddress("tkr_trgseq",tkr_trgseq);
   fChain->SetBranchAddress("cal_trgseq",cal_trgseq);
   fChain->SetBranchAddress("cal_in_mask",cal_in_mask);
   fChain->SetBranchAddress("tkr_layer_enable_0",tkr_layer_enable_0);
   fChain->SetBranchAddress("tkr_layer_enable_1",tkr_layer_enable_1);
   fChain->SetBranchAddress("tkr_out_mask",tkr_out_mask);
   fChain->SetBranchAddress("low",low);
   fChain->SetBranchAddress("high",high);
   fChain->SetBranchAddress("trig_enable",trig_enable);
   fChain->SetBranchAddress("data_mask",data_mask);
   fChain->SetBranchAddress("calib_mask",calib_mask);
   fChain->SetBranchAddress("threshold",threshold);
   fChain->SetBranchAddress("injection",injection);
   fChain->SetBranchAddress("tcc_trg_align",tcc_trg_align);
   fChain->SetBranchAddress("tcc_configuration",tcc_configuration);
   fChain->SetBranchAddress("input_mask",input_mask);
   fChain->SetBranchAddress("trc_csr",trc_csr);
   fChain->SetBranchAddress("ccc_trg_alignment",ccc_trg_alignment);
   fChain->SetBranchAddress("ccc_configuration",ccc_configuration);
   fChain->SetBranchAddress("layer_mask_0",layer_mask_0);
   fChain->SetBranchAddress("layer_mask_1",layer_mask_1);
   fChain->SetBranchAddress("crc_dac",crc_dac);
   fChain->SetBranchAddress("config",config);
   fChain->SetBranchAddress("delay_1",delay_1);
   fChain->SetBranchAddress("delay_2",delay_2);
   fChain->SetBranchAddress("delay_3",delay_3);
   fChain->SetBranchAddress("config_0",config_0);
   fChain->SetBranchAddress("config_1",config_1);
   fChain->SetBranchAddress("ref_dac",ref_dac);
   fChain->SetBranchAddress("log_acpt",log_acpt);
   fChain->SetBranchAddress("fle_dac",fle_dac);
   fChain->SetBranchAddress("fhe_dac",fhe_dac);
   fChain->SetBranchAddress("rng_uld_dac",rng_uld_dac);
   fChain->SetBranchAddress("window_width",&window_width);
   fChain->SetBranchAddress("window_open_mask",&window_open_mask);
   fChain->SetBranchAddress("periodic_prescale",&periodic_prescale);
   fChain->SetBranchAddress("periodic_limit",&periodic_limit);
   fChain->SetBranchAddress("periodic_free",&periodic_free);
   fChain->SetBranchAddress("periodic_onepps",&periodic_onepps);
   fChain->SetBranchAddress("conditions",conditions);
   fChain->SetBranchAddress("engine",engine);
   fChain->SetBranchAddress("towers",towers);
   fChain->SetBranchAddress("acd_cno",&acd_cno);
   fChain->SetBranchAddress("tiles",tiles);
   fChain->SetBranchAddress("tower_busy",&tower_busy);
   fChain->SetBranchAddress("external",&external);
   fChain->SetBranchAddress("configuration",&configuration);
   fChain->SetBranchAddress("r",r);

   fChain->SetBranchAddress("hold_delay_bcast",&hold_delay_bcast);
   fChain->SetBranchAddress("hitmap_delay_bcast",&hitmap_delay_bcast);
   fChain->SetBranchAddress("hitmap_width_bcast",&hitmap_width_bcast);
   fChain->SetBranchAddress("hitmap_deadtime_bcast",&hitmap_deadtime_bcast);
   fChain->SetBranchAddress("veto_delay_bcast",&veto_delay_bcast);
   fChain->SetBranchAddress("veto_width_bcast",&veto_width_bcast);
   fChain->SetBranchAddress("adc_tacq_bcast",&adc_tacq_bcast);
   fChain->SetBranchAddress("pha_threshold_bcast",&pha_threshold_bcast);
   fChain->SetBranchAddress("config_reg_bcast",&config_reg_bcast);
   fChain->SetBranchAddress("tci_dac_bcast",&tci_dac_bcast);
   fChain->SetBranchAddress("bias_dac_bcast",&bias_dac_bcast);
   fChain->SetBranchAddress("hld_dac_bcast",&hld_dac_bcast);
   fChain->SetBranchAddress("veto_dac_bcast",&veto_dac_bcast);
   fChain->SetBranchAddress("veto_vernier_bcast",&veto_vernier_bcast);
   fChain->SetBranchAddress("data_masks_bcast",&data_masks_bcast);
   fChain->SetBranchAddress("tkr_trgseq_bcast",&tkr_trgseq_bcast);
   fChain->SetBranchAddress("cal_trgseq_bcast",&cal_trgseq_bcast);
   fChain->SetBranchAddress("cal_in_mask_bcast",&cal_in_mask_bcast);
   fChain->SetBranchAddress("tkr_layer_enable_0_bcast",&tkr_layer_enable_0_bcast);
   fChain->SetBranchAddress("tkr_layer_enable_1_bcast",&tkr_layer_enable_1_bcast);
   fChain->SetBranchAddress("tkr_out_mask_bcast",&tkr_out_mask_bcast);
   fChain->SetBranchAddress("low_bcast",&low_bcast);
   fChain->SetBranchAddress("high_bcast",&high_bcast);
   fChain->SetBranchAddress("trig_enable_bcast",&trig_enable_bcast);
   fChain->SetBranchAddress("data_mask_bcast",&data_mask_bcast);
   fChain->SetBranchAddress("calib_mask_bcast",&calib_mask_bcast);
   fChain->SetBranchAddress("threshold_bcast",&threshold_bcast);
   fChain->SetBranchAddress("injection_bcast",&injection_bcast);
   fChain->SetBranchAddress("tcc_trg_align_bcast",&tcc_trg_align_bcast);
   fChain->SetBranchAddress("tcc_configuration_bcast",&tcc_configuration_bcast);
   fChain->SetBranchAddress("input_mask_bcast",&input_mask_bcast);
   fChain->SetBranchAddress("trc_csr_bcast",&trc_csr_bcast);
   fChain->SetBranchAddress("ccc_trg_alignment_bcast",&ccc_trg_alignment_bcast);
   fChain->SetBranchAddress("ccc_configuration_bcast",&ccc_configuration_bcast);
   fChain->SetBranchAddress("layer_mask_0_bcast",&layer_mask_0_bcast);
   fChain->SetBranchAddress("layer_mask_1_bcast",&layer_mask_1_bcast);
   fChain->SetBranchAddress("config_bcast",&config_bcast);
   fChain->SetBranchAddress("crc_dac_bcast",&crc_dac_bcast);
   fChain->SetBranchAddress("delay_1_bcast",&delay_1_bcast);
   fChain->SetBranchAddress("delay_2_bcast",&delay_2_bcast);
   fChain->SetBranchAddress("delay_3_bcast",&delay_3_bcast);
   fChain->SetBranchAddress("config_0_bcast",&config_0_bcast);
   fChain->SetBranchAddress("config_1_bcast",&config_1_bcast);
   fChain->SetBranchAddress("ref_dac_bcast",&ref_dac_bcast);
   fChain->SetBranchAddress("log_acpt_bcast",&log_acpt_bcast);
   fChain->SetBranchAddress("fle_dac_bcast",&fle_dac_bcast);
   fChain->SetBranchAddress("fhe_dac_bcast",&fhe_dac_bcast);
   fChain->SetBranchAddress("rng_uld_dac_bcast",&rng_uld_dac_bcast);

   Notify();
}

Bool_t ConfigBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.

   // Get branch pointers
   b_aem_configuration = fChain->GetBranch("aem_configuration");
   b_trgseq = fChain->GetBranch("trgseq");
   b_hold_delay = fChain->GetBranch("hold_delay");
   b_hitmap_delay = fChain->GetBranch("hitmap_delay");
   b_hitmap_width = fChain->GetBranch("hitmap_width");
   b_hitmap_deadtime = fChain->GetBranch("hitmap_deadtime");
   b_veto_delay = fChain->GetBranch("veto_delay");
   b_veto_width = fChain->GetBranch("veto_width");
   b_adc_tacq = fChain->GetBranch("adc_tacq");
   b_pha_threshold = fChain->GetBranch("pha_threshold");
   b_config_reg = fChain->GetBranch("config_reg");
   b_tci_dac = fChain->GetBranch("tci_dac");
   b_bias_dac = fChain->GetBranch("bias_dac");
   b_hld_dac = fChain->GetBranch("hld_dac");
   b_veto_dac = fChain->GetBranch("veto_dac");
   b_veto_vernier = fChain->GetBranch("veto_vernier");
   b_data_masks = fChain->GetBranch("data_masks");
   b_tkr_trgseq = fChain->GetBranch("tkr_trgseq");
   b_cal_trgseq = fChain->GetBranch("cal_trgseq");
   b_cal_in_mask = fChain->GetBranch("cal_in_mask");
   b_tkr_layer_enable_0 = fChain->GetBranch("tkr_layer_enable_0");
   b_tkr_layer_enable_1 = fChain->GetBranch("tkr_layer_enable_1");
   b_tkr_out_mask = fChain->GetBranch("tkr_out_mask");
   b_low = fChain->GetBranch("low");
   b_high = fChain->GetBranch("high");
   b_trig_enable = fChain->GetBranch("trig_enable");
   b_data_mask = fChain->GetBranch("data_mask");
   b_calib_mask = fChain->GetBranch("calib_mask");
   b_threshold = fChain->GetBranch("threshold");
   b_injection = fChain->GetBranch("injection");
   b_tcc_trg_align = fChain->GetBranch("tcc_trg_align");
   b_tcc_configuration = fChain->GetBranch("tcc_configuration");
   b_input_mask = fChain->GetBranch("input_mask");
   b_trc_csr = fChain->GetBranch("trc_csr");
   b_ccc_trg_alignment = fChain->GetBranch("ccc_trg_alignment");
   b_ccc_configuration = fChain->GetBranch("ccc_configuration");
   b_layer_mask_0 = fChain->GetBranch("layer_mask_0");
   b_layer_mask_1 = fChain->GetBranch("layer_mask_1");
   b_crc_dac = fChain->GetBranch("crc_dac");
   b_config = fChain->GetBranch("config");
   b_delay_1 = fChain->GetBranch("delay_1");
   b_delay_2 = fChain->GetBranch("delay_2");
   b_delay_3 = fChain->GetBranch("delay_3");
   b_config_0 = fChain->GetBranch("config_0");
   b_config_1 = fChain->GetBranch("config_1");
   b_ref_dac = fChain->GetBranch("ref_dac");
   b_log_acpt = fChain->GetBranch("log_acpt");
   b_fle_dac = fChain->GetBranch("fle_dac");
   b_fhe_dac = fChain->GetBranch("fhe_dac");
   b_rng_uld_dac = fChain->GetBranch("rng_uld_dac");
   b_window_width = fChain->GetBranch("window_width");
   b_window_open_mask = fChain->GetBranch("window_open_mask");
   b_periodic_prescale = fChain->GetBranch("periodic_prescale");
   b_periodic_limit = fChain->GetBranch("periodic_limit");
   b_periodic_free = fChain->GetBranch("periodic_free");
   b_periodic_onepps = fChain->GetBranch("periodic_onepps");
   b_conditions = fChain->GetBranch("conditions");
   b_engine = fChain->GetBranch("engine");
   b_towers = fChain->GetBranch("towers");
   b_tiles = fChain->GetBranch("tiles");
   b_acd_cno = fChain->GetBranch("acd_cno");
   b_tower_busy = fChain->GetBranch("tower_busy");
   b_external = fChain->GetBranch("external");
   b_configuration = fChain->GetBranch("configuration");
   b_r = fChain->GetBranch("r");

   b_hold_delay_bcast = fChain->GetBranch("hold_delay_bcast");
   b_hitmap_delay_bcast = fChain->GetBranch("hitmap_delay_bcast");
   b_hitmap_width_bcast = fChain->GetBranch("hitmap_width_bcast");
   b_hitmap_deadtime_bcast = fChain->GetBranch("hitmap_deadtime_bcast");
   b_veto_delay_bcast = fChain->GetBranch("veto_delay_bcast");
   b_veto_width_bcast = fChain->GetBranch("veto_width_bcast");
   b_adc_tacq_bcast = fChain->GetBranch("adc_tacq_bcast");
   b_pha_threshold_bcast = fChain->GetBranch("pha_threshold_bcast");
   b_config_reg_bcast = fChain->GetBranch("config_reg_bcast");
   b_tci_dac_bcast = fChain->GetBranch("tci_dac_bcast");
   b_bias_dac_bcast = fChain->GetBranch("bias_dac_bcast");
   b_hld_dac_bcast = fChain->GetBranch("hld_dac_bcast");
   b_veto_dac_bcast = fChain->GetBranch("veto_dac_bcast");
   b_veto_vernier_bcast = fChain->GetBranch("veto_vernier_bcast");
   b_data_masks_bcast = fChain->GetBranch("data_masks_bcast");
   b_tkr_trgseq_bcast = fChain->GetBranch("tkr_trgseq_bcast");
   b_cal_trgseq_bcast = fChain->GetBranch("cal_trgseq_bcast");
   b_cal_in_mask_bcast = fChain->GetBranch("cal_in_mask_bcast");
   b_tkr_layer_enable_0_bcast = fChain->GetBranch("tkr_layer_enable_0_bcast");
   b_tkr_layer_enable_1_bcast = fChain->GetBranch("tkr_layer_enable_1_bcast");
   b_tkr_out_mask_bcast = fChain->GetBranch("tkr_out_mask_bcast");
   b_low_bcast = fChain->GetBranch("low_bcast");
   b_high_bcast = fChain->GetBranch("high_bcast");
   b_trig_enable_bcast = fChain->GetBranch("trig_enable_bcast");
   b_data_mask_bcast = fChain->GetBranch("data_mask_bcast");
   b_calib_mask_bcast = fChain->GetBranch("calib_mask_bcast");
   b_threshold_bcast = fChain->GetBranch("threshold_bcast");
   b_injection_bcast = fChain->GetBranch("injection_bcast");
   b_tcc_trg_align_bcast = fChain->GetBranch("tcc_trg_align_bcast");
   b_tcc_configuration_bcast = fChain->GetBranch("tcc_configuration_bcast");
   b_input_mask_bcast = fChain->GetBranch("input_mask_bcast");
   b_trc_csr_bcast = fChain->GetBranch("trc_csr_bcast");
   b_ccc_trg_alignment_bcast = fChain->GetBranch("ccc_trg_alignment_bcast");
   b_ccc_configuration_bcast = fChain->GetBranch("ccc_configuration_bcast");
   b_layer_mask_0_bcast = fChain->GetBranch("layer_mask_0_bcast");
   b_layer_mask_1_bcast = fChain->GetBranch("layer_mask_1_bcast");
   b_crc_dac_bcast = fChain->GetBranch("crc_dac_bcast");
   b_config_bcast = fChain->GetBranch("config_bcast");
   b_delay_1_bcast = fChain->GetBranch("delay_1_bcast");
   b_delay_2_bcast = fChain->GetBranch("delay_2_bcast");
   b_delay_3_bcast = fChain->GetBranch("delay_3_bcast");
   b_config_0_bcast = fChain->GetBranch("config_0_bcast");
   b_config_1_bcast = fChain->GetBranch("config_1_bcast");
   b_ref_dac_bcast = fChain->GetBranch("ref_dac_bcast");
   b_log_acpt_bcast = fChain->GetBranch("log_acpt_bcast");
   b_fle_dac_bcast = fChain->GetBranch("fle_dac_bcast");
   b_fhe_dac_bcast = fChain->GetBranch("fhe_dac_bcast");
   b_rng_uld_dac_bcast = fChain->GetBranch("rng_uld_dac_bcast");

   return kTRUE;
}

void ConfigBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ConfigBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef ConfigBase_cxx

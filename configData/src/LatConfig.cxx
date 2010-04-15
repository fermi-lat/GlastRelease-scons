
#include "configData/LatConfig.h"

#include "TTree.h"
#include <iostream>


TTree* LatConfig::makeTree(const std::string& prefix) {
  if ( m_tree != 0 ) {
    delete m_tree;
  }
  m_tree = new TTree(m_config.name().c_str(),"Configuration Tree");
  m_config.makeBranches(*m_tree,prefix);
  m_gem.makeBranches(*m_tree,prefix);
  return m_tree;
}

void LatConfig::latch() {
  if ( m_tree == 0 ) return;
  m_tree->Fill();
}

void LatConfig::initConfig() {
  
  // AFE
  ChannelKey afeSize(12,18,1,1);
  ChannelKey arcSize(12,1,1,1);
  ChannelKey singleton(1,1,1,1);

  ChannelKey temSize(16,1,1,1);
  ChannelKey sptSize(16,36,1,1);
  ChannelKey tfeSize(16,36,24,1);
  ChannelKey tccSize(16,8,1,1);
  ChannelKey trcSize(16,8,9,1);

  ChannelKey cccSize(16,4,1,1);
  ChannelKey crcSize(16,4,4,1);
  ChannelKey cfeSize(16,4,4,12);

  // AEM
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("aem_configuration",'i',singleton)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("trgseq",'i',singleton)));

  // ARC
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("hold_delay",'s',arcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("hitmap_delay",'s',arcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("hitmap_width",'s',arcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("hitmap_deadtime",'s',arcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("veto_delay",'s',arcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("veto_width",'s',arcSize)));    
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("adc_tacq",'s',arcSize)));    
  // ARC Special
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("pha_threshold",'s',afeSize)));


  // AFE
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("config_reg",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("tci_dac",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("bias_dac",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("hld_dac",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("veto_dac",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("veto_vernier",'s',afeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UShort_t>("acd_id",'s',afeSize)));

  // TEM
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("data_masks",'i',temSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tkr_trgseq",'i',temSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("cal_trgseq",'i',temSize)));
  
  // TIC
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("cal_in_mask",'i',temSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tkr_layer_enable_0",'i',temSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tkr_layer_enable_1",'i',temSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tkr_out_mask",'i',temSize)));
  

  // SPT
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("low",'i',sptSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("high",'i',sptSize)));

  // TFE
  m_config.addBranch( *(new ConfigBranchImpl<ULong64_t>("trig_enable",'l',tfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<ULong64_t>("data_mask",'l',tfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<ULong64_t>("calib_mask",'l',tfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("threshold",'i',tfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("injection",'i',tfeSize)));

  // TCC
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tcc_trg_align",'i',tccSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("tcc_configuration",'i',tccSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("input_mask",'i',tccSize)));

  // TRC
  m_config.addBranch( *(new ConfigBranchImpl<ULong64_t>("trc_csr",'l',trcSize)));

  // CCC
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("ccc_trg_alignment",'i',cccSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("ccc_configuration",'i',cccSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("layer_mask_0",'i',cccSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("layer_mask_1",'i',cccSize)));

  // CRC
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("config",'i',crcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("crc_dac",'i',crcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("delay_1",'i',crcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("delay_2",'i',crcSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("delay_3",'i',crcSize)));
 
  // CFE
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("config_0",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("config_1",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("ref_dac",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("log_acpt",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("fle_dac",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("fhe_dac",'i',cfeSize)));
  m_config.addBranch( *(new ConfigBranchImpl<UInt_t>("rng_uld_dac",'i',cfeSize)));
   
}

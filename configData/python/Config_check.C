#define Config_check_cxx
#include "Config_check.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iomanip>

void Config_check::TKR_Timing()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   assert(nentries==1);
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   startcolor(black);
   std::cout<<"File: "<<m_filename<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<"+             Config                 +"<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<std::endl;
   std::cout<<"TEM tcc_trg_align"<<std::endl;
   std::cout<<"----------------------"<<std::endl;
   switchcolor(magenta);
   std::cout<<std::setw(3) <<"BC"<<hex<<setw(10)<<tcc_trg_align_bcast<<dec<<std::endl;
   for (int i=0;i<16;i++){
     for (int j=1;j<8;j++) assert(tcc_trg_align[i][j]==tcc_trg_align[i][0]);
     switchcolor(black);
     std::cout<<std::setw(3)<<i;
     if(tcc_trg_align[i][0]==tcc_trg_align_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<hex<<setw(10)<<tcc_trg_align[i][0]<<dec<<std::endl;
   }
   switchcolor(black);
   std::cout<<std::endl;
   std::cout<<"It has been verified that all tcc_trg_align registers in a tower are identical."<<std::endl;
      // if (Cut(ientry) < 0) continue;
   std::cout<<std::endl;
   std::cout<<"Legend: "<<std::endl;
   std::cout<<color(magenta)<<"Magenta: Broadcast value"<<std::endl;
   endcolor();
   std::cout<<color(green)<<"Green: Register value that is same as broadcast"<<std::endl;
   endcolor();
   std::cout<<color(red)<<"Red: Register value that is different from broadcast"<<color(black)<<std::endl;
   endcolor();
   std::cout<<std::endl;
}
void Config_check::ACD_Timing()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   assert(nentries==1);
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   startcolor(black);
   std::cout<<"File: "<<m_filename<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<"+             ACD_Timing                 +"<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<std::endl;
   std::cout<<"GARC veto_delay veto_width hitmap_delay hitmap_width hitmap_deadtime hold_delay adc_tacq"<<std::endl;
   std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
   switchcolor(magenta);
   std::cout<<std::setw(3)<<"BC"<<setw(8)<<veto_delay_bcast<<setw(10)<<veto_width_bcast;
   std::cout<<setw(12)<<hitmap_delay_bcast<<setw(12)<<hitmap_width_bcast;
   std::cout<<setw(13)<<hitmap_deadtime_bcast<<setw(16)<<hold_delay_bcast;
   std::cout<<setw(11)<<adc_tacq_bcast<<std::endl;
   for (int i=0;i<12;i++){
     switchcolor(black);
     std::cout<<std::setw(3)<<i;
     if(veto_delay[i]==veto_delay_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(8)<<veto_delay[i];
     if(veto_width[i]==veto_width_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(10)<<veto_width[i];
     if(hitmap_delay[i]==hitmap_delay_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(12)<<hitmap_delay[i];
     if(hitmap_width[i]==hitmap_width_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(12)<<hitmap_width[i];
     if(hitmap_deadtime[i]==hitmap_deadtime_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(13)<<hitmap_deadtime[i];
     if(hold_delay[i]==hold_delay_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(16)<<hold_delay[i];
     if(adc_tacq[i]==adc_tacq_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(11)<<adc_tacq[i]<<std::endl;
}
   switchcolor(black);
   std::cout<<std::endl;
      // if (Cut(ientry) < 0) continue;
   std::cout<<std::endl;
   std::cout<<"Legend: "<<std::endl;
   std::cout<<color(magenta)<<"Magenta: Broadcast value"<<std::endl;
   endcolor();
   std::cout<<color(green)<<"Green: Register value that is same as broadcast"<<std::endl;
   endcolor();
   std::cout<<color(red)<<"Red: Register value that is different from broadcast"<<color(black)<<std::endl;
   endcolor();
   std::cout<<std::endl;
}
void Config_check::CAL_Timing()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   assert(nentries==1);
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   startcolor(black);
   std::cout<<"File: "<<m_filename<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<"+             CAL_Timing                 +"<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<std::endl;
   std::cout<<"TEM ccc_trg_alignment delay_1 delay_2 delay_3"<<std::endl;
   std::cout<<"---------------------------------------------------"<<std::endl;
   switchcolor(magenta);
   std::cout<<std::setw(3)<<"BC"<<hex<<setw(15)<<ccc_trg_alignment_bcast<<dec<<setw(10)<<delay_1_bcast<<setw(8)<<delay_2_bcast<<setw(8)<<delay_3_bcast<<std::endl;
   for (int i=0;i<16;i++){
     for (int j=0;j<4;j++){
       for (int k=1;k<4;k++){
	 assert (delay_1[i][j][k]==delay_1[i][0][0]);
	 assert (delay_2[i][j][k]==delay_2[i][0][0]);
	 assert (delay_3[i][j][k]==delay_3[i][0][0]);
       }
       assert (ccc_trg_alignment[i][j]==ccc_trg_alignment[i][0]);
     }
     switchcolor(black);
     std::cout<<std::setw(3)<<i;
     if(ccc_trg_alignment[i][0]==ccc_trg_alignment_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<hex<<setw(15)<<ccc_trg_alignment[i][0];
     if(delay_1[i][0][0]==delay_1_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<dec<<setw(10)<<delay_1[i][0][0];
     if(delay_2[i][0][0]==delay_2_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(8)<<delay_2[i][0][0];
     if(delay_3[i][0][0]==delay_3_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<setw(8)<<delay_3[i][0][0]<<std::endl;
   }
   switchcolor(black);
   std::cout<<std::endl;
   std::cout<<"All registers of a type have the identical value within one TEM"<<std::endl;
      // if (Cut(ientry) < 0) continue;
   std::cout<<std::endl;
   std::cout<<"Legend: "<<std::endl;
   std::cout<<color(magenta)<<"Magenta: Broadcast value"<<std::endl;
   endcolor();
   std::cout<<color(green)<<"Green: Register value that is same as broadcast"<<std::endl;
   endcolor();
   std::cout<<color(red)<<"Red: Register value that is different from broadcast"<<color(black)<<std::endl;
   endcolor();
   std::cout<<std::endl;
}
void Config_check::GNL_Timing()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   assert(nentries==1);
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   startcolor(black);
   std::cout<<"File: "<<m_filename<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<"+             GNL_Timing                 +"<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<std::endl;
   std::cout<<"AEM trgseq: "<<trgseq<<std::endl;
   std::cout<<std::endl;
   std::cout<<"TEM cal_trgseq tkr_trgseq"<<std::endl;
   std::cout<<"-------------------------"<<std::endl;
   color col(magenta);
   std::cout<<col;
   std::cout<<std::setw(3)<<"BC"<<setw(8)<<cal_trgseq_bcast<<setw(11)<<tkr_trgseq_bcast<<std::endl;
   col.fore=black;
   std::cout<<col;
   for (int i=0;i<16;i++){
     std::cout<<std::setw(3)<<i;
     if (cal_trgseq[i]==cal_trgseq_bcast)col.fore=green;
     else col.fore=red;
     std::cout<<col;
     std::cout<<setw(8)<<cal_trgseq[i];
     if (tkr_trgseq[i]==tkr_trgseq_bcast)col.fore=green;
     else col.fore=red;
     std::cout<<col;
     std::cout<<setw(11)<<tkr_trgseq[i];
     col.fore=black;
     std::cout<<col;
     std::cout<<std::endl;
   }
   std::cout<<std::endl;
   std::cout<<"Legend: "<<std::endl;
   std::cout<<color(magenta)<<"Magenta: Broadcast value"<<std::endl;
   endcolor();
   std::cout<<color(green)<<"Green: Register value that is same as broadcast"<<std::endl;
   endcolor();
   std::cout<<color(red)<<"Red: Register value that is different from broadcast"<<color(black)<<std::endl;
   endcolor();
   std::cout<<std::endl;
      // if (Cut(ientry) < 0) continue;
}

void Config_check::GNL_Mode()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   assert(nentries==1);
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   startcolor(black);
   std::cout<<"File: "<<m_filename<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<"+             GNL_Mode                   +"<<std::endl;
   std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
   std::cout<<std::endl;
   std::cout<<"AEM configuration: "<<hex<<aem_configuration<<dec<<std::endl;
   std::cout<<std::endl;
   std::cout<<"TEM   data_masks cal_in_mask tkr_out_mask tkr_layer_enable_0 tkr_layer_enable_1"<<std::endl;
   std::cout<<"-------------------------------------------------------------------------------"<<std::endl;
   color col;
   col.fore=magenta;
   std::cout<<col; 
   std::cout<<std::setw(4)<<"BC"<<hex<<setw(10)<<data_masks_bcast<<setw(12)<<cal_in_mask_bcast;
   std::cout<<std::setw(14)<<tkr_out_mask_bcast;
   std::cout<<std::setw(16)<<tkr_layer_enable_0_bcast;
   std::cout<<std::setw(18)<<tkr_layer_enable_1_bcast;
   std::cout<<endl;
   switchcolor(black);
   for (int i=0;i<16;i++){
     
     std::cout<<std::setw(4)<<i;
     if(data_masks[i]==data_masks_bcast)switchcolor(green);
     else switchcolor(red);
     std::cout<<hex<<setw(10)<<data_masks[i];
     if(cal_in_mask[i]==cal_in_mask_bcast)switchcolor(green);
     else col.fore=red;
     std::cout<<hex<<setw(12)<<cal_in_mask[i];
     if(tkr_out_mask[i]==tkr_out_mask_bcast)switchcolor(green);
     else col.fore=red;
     std::cout<<std::setw(14)<<tkr_out_mask[i];
     if(tkr_layer_enable_0[i]==tkr_layer_enable_0_bcast)switchcolor(green);
     else col.fore=red;
     std::cout<<std::setw(16)<<tkr_layer_enable_0[i];
     if(tkr_layer_enable_1[i]==tkr_layer_enable_1_bcast)switchcolor(green);
     else col.fore=red;
     std::cout<<std::setw(18)<<tkr_layer_enable_1[i];
     switchcolor(black);
     std::cout<<dec<<std::endl;
   }
   std::cout<<std::endl;
   std::cout<<"Legend: "<<std::endl;
   col.fore=magenta;
   std::cout<<col<<"Magenta: Broadcast value"<<std::endl;
   endcolor();
   col.fore=green;
   std::cout<<col<<"Green: Register value that is same as broadcast"<<std::endl;
   endcolor();
   col.fore=red;
   std::cout<<col<<"Red: Register value that is different from broadcast"<<color(black)<<std::endl;
   col.fore=black;
   std::cout<<col;
   endcolor();
   std::cout<<std::endl;
}

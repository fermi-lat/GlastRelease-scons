/*! Trigger Analysis Histogram Defintions */

void RootTreeAnalysis::TriggerHistDefine() {
    
    histFile->cd();

    TH1F *TRGBITS = new TH1F("TRGBITS", "Trigger Bits Set",
        10, 0, 10);
    TRGBITS->SetXTitle("Bit");

    TH1F *TRGWORD = new TH1F("TRGWORD", "Trigger Word",
        32, 0, 32);
    TRGWORD->SetXTitle("Word Value");
    
}

/* Process the Monte Carlo Data
Called by Go()
*/


void RootTreeAnalysis::Trigger() {
    
// get multiplicities

     L1T l1T = evt->getL1T();

     int trigWord = l1T->getTriggerWord(); 
     ((TH1F*)GetObjectPtr("TRGWORD"))->Fill(trigWord);

     for (int ib=0; ib<L1T::number_of_trigger_bits; ib++) {
	 if (trigWord & 1<<ib) ((TH1F*)GetObjectPtr("TRGBITS"))->Fill(ib);
     }

}


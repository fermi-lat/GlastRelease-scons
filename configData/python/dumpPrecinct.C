void dumpPrecinct(const char* filename, const char* myMacro, const char* precinct) 
{
  gROOT->LoadMacro(myMacro);
  Config_check f(filename);
  cout << "Precinct: " << precinct << endl;
  if (precinct == "GNL_Mode")
    f.GNL_Mode();
  else if (precinct == "GNL_Timing")
    f.GNL_Timing();
  else if (precinct == "TKR_Timing")
    f.TKR_Timing();
  else if (precinct == "ACD_Timing")
    f.ACD_Timing();
  else if (precinct == "CAL_Timing")
    f.CAL_Timing();
  else
    cout << "couldn't understand " << precinct << endl;
}

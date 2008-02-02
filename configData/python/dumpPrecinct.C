void dumpPrecinct(const char* filename, const char* precinct) 
{
  cout << precinct << endl;
  gROOT->LoadMacro("Config_check.C");
  Config_check f(filename);
  switch (*precinct) {
  case 'GNL_Timing':
    f.GNL_Timing();
    break;
  case 'TKR_Timing':
    f.TKR_Timing();
    break;
  case 'ACD_Timing':
    f.ACD_Timing();
    break;
  case 'CAL_Timing':
    f.CAL_Timing();
    break;
  default:
    cout << "couldn't understand " << precinct << endl;
  }
}

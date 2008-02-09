void dumpPrecinct(const char* filename, const char* myMacro, const char* precinct) 
{
  gROOT->LoadMacro(myMacro);
  Config_check f(filename);
  switch (*precinct) {
  case 'GNL_Mode':
    f.GNL_Mode();
    break;
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

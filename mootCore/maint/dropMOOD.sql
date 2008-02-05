# Drop all MOOD tables in case we want to start over

####use mood_test;

## drop relation tables first.  Comment out drop commands for planned tables
## that do not exist yet

drop table if exists Configs_to_FSW;
drop table if exists FSW_to_Parameters;
drop table if exists Configs_to_Parameters;

drop table if exists FClass_to_PClass;
drop table if exists OClass_to_AClass;
drop table if exists Vote_PClass_AClass;
drop table if exists Parameters_to_Ancillary;
drop table if exists Offline_to_Ancillary;
drop table if exists Container_Precinct;

### drop alias tables
drop table if exists Vote_aliases;
drop table if exists Ancillary_aliases;

### drop grouping tables
drop table if exists Acq_mode;
drop table if exists Acq_summary;
drop table if exists Configs;

### drop instance tables
drop table if exists FSW_inputs;
drop table if exists Parameters;
drop table if exists Votes;
drop table if exists Ancillary;
drop table if exists OfflineCalib;

drop table if exists LATC_Backup;
drop table if exists Parameter_class;

drop table if exists Precincts;
drop table if exists Instrument;

###drop class tables
drop table if exists FSW_class;
drop table if exists Offline_class;
drop table if exists Ancillary_class;

drop table if exists Delegates;


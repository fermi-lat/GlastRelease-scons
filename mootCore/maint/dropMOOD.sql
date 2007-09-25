# Drop all MOOD tables in case we want to start over

####use mood_test;

## drop relation tables first.  Comment out drop commands for planned tables
## that do not exist yet

drop table Configs_to_FSW;
drop table FSW_to_Parameters;
drop table Configs_to_Parameters;

drop table FClass_to_PClass;
drop table OClass_to_AClass;
drop table Vote_PClass_AClass;
drop table Parameters_to_Ancillary;
drop table Offline_to_Ancillary;
drop table Container_Precinct;

### drop grouping tables
drop table Configs;

### drop alias tables
drop table Vote_aliases;
drop table Ancillary_aliases;



### drop instance tables
drop table FSW_inputs;
drop table Parameters;
drop table Votes;
drop table Ancillary;
drop table OfflineCalib;

drop table LATC_Backup;
drop table Parameter_class;

drop table Precincts;
drop table Instrument;

###drop class tables
drop table FSW_class;
drop table Offline_class;
drop table Ancillary_class;

drop table Delegates;

# Drop all MOOD tables in case we want to start over

####use mood_test;

## drop relation tables first.  Comment out drop commands for planned tables
## that don't exist yet

###drop table History_to_Parameters.sql;
drop table Configs_to_FSW;
drop table FSW_to_Parameters;
### drop table Master_to_LATC;  probably don't need this one
drop table FClass_to_PClass;
drop table PClass_to_AClass;
drop table OClass_to_AClass;
drop table Parameters_to_Ancillary;
drop table Offline_to_Ancillary;
##drop table History_to_Parameters;


###drop table History;

### drop grouping tables
drop table Configs;
## drop table CalibSet;

### drop instance tables
drop table FSW_inputs;
drop table Parameters;
drop table Ancillary;
drop table OfflineCalib;

###drop class tables
drop table FSW_class;
drop table Parameter_class;
drop table Offline_class;
drop table Ancillary_class;


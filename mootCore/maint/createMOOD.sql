# Create tables for MOOD.  Order so that primary keys are defined
# before the foreign keys that refer to them.

# for the time being, use only test dbs

## use mood_test;
## assume database has already been selected.  For example, can define
## environment variable for dbname (say thedb), and write a 2-line
## script to a. do the connection using ${thedb} and b. invoke this
## sql file.

##  Seems to be no place for a history table; better to get such information
## from start run commands
##  source createHistory.sql;

## This file creates all the tables needed for a particular MOOD-style
## database.  There will probably be one per site of interest:
##   mood_SLAC
##   mood_NRL
##   mood_MOR   (mobile rack)

## Tables in MOOD are of several types: 

##  Type                Example            Typical row
##  instance table      Parameters         reference to a file
##  class table         Parameter_class    describes class of Parameters
##  relation table      FSW_to_Parameters  relates entries from 2 other tables
##  collection          Config             provide a handle to group together
##                                         instances from another table (here
##                                         FSW_Inputs)
##  history             History            provide handle for time interval
##                                         over which config was fixed

## Some care must be taken about the order in which tables are created or
## dropped since some have foreign keys referring to others. The following
## is an acceptable ordering for creation
##    Type                                                     Count
## 1. Class tables                                             4/4
## 2. Instance tables (usually refer to a class table)         4/4
## 3. Collection tables                                        1/2
## 4. History                                                  0/1
## 5. Relation tables (may refer to any of the other types)    7/9
## In m/n n=number planned, m is number fully specified
##
## Tables may be dropped in reverse order

## Delegates must be created before Parameter_class and Conditions
## Both will some day refer to it.
source createDelegates.sql;      

## class tables
source createFSW_class.sql;
source createParameter_class.sql;
source createOffline_class.sql;   
source createAncillary_class.sql; 

## instance tables
## Conditions must be created before Parameters since the latter
## will have (some day) a foreign key referencing Conditions.
source createConditions.sql;

source createParameters.sql;
source createFSW_inputs.sql;
source createAncillary.sql;   
source createOfflineCalib.sql;  

## collections
source createConfigs.sql;


# Supporting tables; mostly relation tables
source createFClass_to_PClass.sql;
source createConfigs_to_FSW.sql;
source createFSW_to_Parameters.sql;
source createConfigs_to_Parameters.sql;

source createParameters_to_Ancillary.sql;   
source createOffline_to_Ancillary.sql;

source createPClass_to_AClass.sql; 
source createOClass_to_AClass.sql; 

## sui generis
source createInstrument.sql;



# $Header$
# File for creating metadata table for calibration
# This (or any other) script may be run from within mysql 
#     mysql>  source createOfflineCalib.sql
# or directly from shell prompt (known as "batch mode"), e.g.
#     $ mysql -p -h [hostname] [dbname] < createOfflineCalib.sql
#       input_desc       Comments about input to this calibration procedure
#       notes            A spot for miscellaneous comments
create table OfflineCalib
  (ser_no int not null auto_increment primary key, 
   instrument  varchar(16) not null
      COMMENT 'e.g. LAT, EM2  ',
   calib_type  varchar(20) not null
     COMMENT 'should be one of the types in Offline_class',
   flavor      varchar(30) not null default 'vanilla',
   data_fmt    varchar(10) not null,
   data_ident     varchar(255) not null,
   vstart      datetime not null  default '2037-12-31 00:00:00'
COMMENT 'start of validity interval; how set depends on value of config in class table',
   vend        datetime not null default '2037-12-31 00:00:00',
   proc_level ENUM('PROD', 'TEST', 'DEV', 'SUPSED') not null default 'TEST',
   status                 ENUM('STARTED', 'ABORT', 'CREATED', 'INVALID')
                            NOT NULL DEFAULT 'STARTED'
COMMENT 'refers to status of entry, in case file has to be copied',
   completion ENUM('OK','INC', 'ABORT') not null,
   locale      varchar(30) not null default 'orbit',
   input_desc  varchar(255),
   notes varchar(255),
   fmt_version varchar(16),
   update_time  datetime,
   enter_time  datetime,
   prod_start  datetime,
   prod_end    datetime,
   creator     varchar(255),
   uid         varchar(12) not null,
   data_size   int(11),
   index typical_search (completion, instrument, calib_type, flavor, proc_level)  ) TYPE=InnoDB   COMMENT="calibration metadata table after move to MOOT db";

#    Following won't be implemented until and unless it's necessary
#   read_status      In case we want to, e.g., lock out readers while a 
#                    row is being created or updated.
###   Probably don't need a 'busy bit' per row
###    read_status ENUM('BUSY', 'OK', 'BAD') not null default 'BUSY', 

# $Header$
# File for creating metadata table for calibration
# This (or any other) script may be run from within mysql 
#     mysql>  source createmeta.sql
# or directly from shell prompt (known as "batch mode"), e.g.
#     $ mysql -p < createmeta.sql


use calib;


## start out with alternate name so as not to disturb current working
## version                    
#    Field meanings:
#       ser_no           serial number.  One series is used for all
#                        calibration types. (automatic assignment)
#       instrument       flight, EM, etc.
#       calib_type       What kind of data, e.g. TKR bad strips
#       flavor           Denotes special kind of calibration; defaults to "vanilla"
#       data_fmt         Anticipate supporting XML and ROOT (just XML to start)
#       data_size        Optional field; interpretation to depend on calib_type
#       vstart,vend      Validity interval for data
#       enter_time       When this record was entered in dbs (automatic assignment)
#       fmt_version      Versioning field placeholder so that code has option
#                        of, e.g., reading in old data.
#       proc_level       How seriously are we to take this data?
#       completion       Did it?     
#       read_status      In case we want to, e.g., lock out readers while a 
#                        row is being created or updated.
#       creator          What program created the data?
#       uid              What person created/entered the data?
#       data_ident       Filename or other identifier for persistent data
# $Header$
#       input_desc       Comments about input to this calibration procedure
#       notes            A spot for miscellaneous comments
create table metadata_v2
  (ser_no mediumint not null auto_increment primary key, 
   instrument  varchar(16) not null,
   calib_type  varchar(20) not null,
   flavor      varchar(30) not null default 'VANILLA',
   data_fmt    varchar(10) not null,
   data_size   int(11),
   vstart      datetime,
   vend        datetime,
   enter_time  timestamp,
   fmt_version varchar(12),
   completion ENUM('OK','INC', 'ABORT'),
   proc_level ENUM('PROD', 'TEST', 'DEV', 'SUPSED') not null default 'TEST',
   creator     varchar(255),
   uid         varchar(12) not null,
   data_ident     varchar(255) not null,  
   input_desc  varchar(255),
   notes varchar(255),
   index typical_search (completion, instrument, calib_type, proc_level)  )
  comment = "calibration metadata table (test version)";

###   Probably don't need a 'busy bit' per row
###    read_status ENUM('BUSY', 'OK', 'BAD') not null default 'BUSY', 

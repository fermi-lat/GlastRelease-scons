#  $Header$
#  File for creating history table.  A row is added to this table every
#  time a structural change is made to the calibration metadata table
#  or a new version of the metadata table is created

use calib;

create table history
 (ser_no mediumint not null auto_increment primary key,
  table_name      varchar(30) not null,
  operation       varchar(30) not null,
  operation_time  datetime,
  enter_time      timestamp,
  details         varchar(255)
 )
 comment = "Log structural changes to calibration metadata table(s?)";

#   ser_no         -- unique label per entry (set automatically)
#   table_name     -- name of calibration table created or modified
#   operation      -- e.g., "create", "add columns", etc.
#   operation_time -- approximate time operation occurred
#   enter_time     -- time this entry was made (value set automatically)
#   details        -- if more info than will fit in operation field 

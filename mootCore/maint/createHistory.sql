create table History
 (history_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  ground_creation          TIMESTAMP
     COMMENT 'creation time of row',
  ground_update            TIMESTAMP
     COMMENT 'last time the same config was seen in data stream',
  flight_start             TIMESTAMP
     COMMENT 'flight timestamp from first data packet with this config',
  flight_update            TIMESTAMP   
     COMMENT 'flight timestamp from most recent data packet with this config',
  flight_interval          FLOAT  NOT NULL DEFAULT 0.0
     COMMENT 'seconds between old flight_update and new',   
  application            VARCHAR(16)      NOT NULL
     COMMENT 'typical value would be FSW, maybe HV.',
  state                  VARCHAR(16)      NOT NULL 
     COMMENT 'One of, e.g., FROM_DATA, REQUEST, REQUEST_ACK; could be enum'
  ids                    VARCHAR(50)
     COMMENT 'list of ids.  Or do we want individual columns for a known set?',
  mode                   VARCHAR(16)
     COMMENT 'not sure if we want to keep track of that here..'              
  INDEX(flight_start)
 ) TYPE=InnoDB
COMMENT='a row corresponds to a collection of logical ids as reported in 
  the datastream. Make a new row when collection changes or time gap is large';


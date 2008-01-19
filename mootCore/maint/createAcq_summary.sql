## Contains per-acquisition configuration information.  See also Acq_mode table, which
## contains information concerning config info which potentially varies on a per-mode basis

create table Acq_summary
 (acq_summary_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
    COMMENT "primary key",
  scid INT UNSIGNED NOT NULL
    COMMENT 'Same as SCID in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY',
  started_at INT UNSIGNED NOT NULL 
    COMMENT 'Same as STARTEDAT in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY',
  config_fk INT UNSIGNED
    COMMENT 'Refers to MOOT config presumed used for this acquisition',
  type VARCHAR(16) NOT NULL
    COMMENT 'Same as TYPE in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY, typically LCI or LPA',
  analysis VARCHAR(64)
    COMMENT 'Same as ANALYSIS in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY; should be derivable
             from config_fk',
  status VARCHAR(16) NOT NULL DEFAULT 'STARTED'
    COMMENT 'Whether acquisition is done, in progress, aborted, etc.',
  
  nevts INT UNSIGNED  NOT NULL DEFAULT '0' 
    COMMENT 'Number of events from this acq decoded so far',
  evtutc0 DATETIME  NOT NULL 
    COMMENT 'timestamp for first decoded event so far for this acq, truncated to nearest second',
  evtutc1 DATETIME  NOT NULL 
    COMMENT 'timestamp for last decoded event so far for this acq, rounded up to nearest second',
  hw_key INT UNSIGNED NOT NULL DEFAULT '0'
    COMMENT 'Same as HWKEY (fmx key for LATC master) in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY',
  sw_key INT UNSIGNED NOT NULL DEFAULT '0'
    COMMENT 'Same as SWKEY (sbs key or lci script key) in Oracle table ISOC_FLIGHT.GLASTOPS_ACQSUMMARY',
  creator                VARCHAR(32) NOT NULL
    COMMENT 'uid of last updater of this record',
  creation_time  DATETIME
    COMMENT 'utc time of initial insert of this row',
  update_time  DATETIME
    COMMENT 'utc time of last update of this row',
  FOREIGN KEY(config_fk) REFERENCES Configs (config_key)
  ON DELETE SET NULL
  ON UPDATE CASCADE,
  UNIQUE INDEX (scid, started_at)
 ) TYPE=InnoDB 
COMMENT='A row corresponds to a commanded acquisition';

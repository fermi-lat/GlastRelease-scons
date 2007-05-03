create table Configs
 (config_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  instrument  varchar(16) not null default 'LAT'
      COMMENT 'e.g. LAT, EM2  ',
  creation_request_time  DATETIME
    COMMENT 'time of initial insert of this row',
  creation_end_time      DATETIME DEFAULT '2037-12-31'
    COMMENT 'time this config was completed or attempt was abandoned',
  status                 ENUM('STARTED', 'ABORT', 'CREATED', 'INVALID')
                            NOT NULL DEFAULT 'STARTED',
  active_state           ENUM('ACTIVE', 'INACTIVE', 'MULTI') NOT NULL
                           DEFAULT 'INACTIVE' 
  COMMENT 'App  preference for Config of particular name, alg,step,instr',
  vote_fk INT UNSIGNED 
    COMMENT 'refers to container Vote file describing intent (optional)',

  creator                VARCHAR(32) NOT NULL,
  mode                   ENUM('DATA', 'CALIB', 'QUIET', 'TEST') NOT NULL
    COMMENT 'determines which parameters will belong to this Config',
  name                   VARCHAR(255) NOT NULL,
  algorithm              VARCHAR(255) NOT NULL DEFAULT 'none'
   COMMENT 'can take several runs, each with different config, to do one calibration', 
  alg_step                  SMALLINT NOT NULL DEFAULT 1 
    COMMENT 'which step if part of multi-step algorithm', 
  description            VARCHAR(255) NOT NULL,
  FOREIGN KEY(vote_fk) REFERENCES Votes (vote_key)
 ) TYPE=InnoDB 
COMMENT='A row corresponds to a particular LAT configuration';


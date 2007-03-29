create table Parameters
 (parm_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,

  class_fk               INT UNSIGNED NOT NULL 
     COMMENT 'refers to row in Parm_class',

  vote_fk               INT UNSIGNED NOT NULL DEFAULT '1'
     COMMENT 'refers to row in Votes',

  instrument  varchar(16) not null default 'LAT'
     COMMENT 'e.g. LAT, EM2  ',

  source                 VARCHAR(255) NOT NULL 
     COMMENT 'path to file in Moot archive, relative to root',

  source_fmt             VARCHAR(32) NOT NULL
     COMMENT 'could be, e.g. xml',
  precursor             VARCHAR(255) NOT NULL DEFAULT 'unknown'
     COMMENT 'source for the source',
  creation_time         DATETIME 
  COMMENT 'When this row was created',

  src_creation_time         DATETIME 
  COMMENT 'When the file this row references  was created',


  status                 ENUM('STARTED', 'ABORT', 'CREATED', 'INVALID')
                            NOT NULL DEFAULT 'STARTED'
COMMENT 'refers to status of entry, in case file has to be copied',

  creator               VARCHAR(32) NOT NULL,

  quality               ENUM('PROD', 'DEV', 'TEST', 'SUPSED', 'INVALID') 
                       NOT NULL,
  description           VARCHAR(255) NOT NULL
     COMMENT 'a place to mention anything special about this parameter',
  checksum              INT UNSIGNED  
     COMMENT 'provides quick way for mostly-reliable diff',
  size                  INT UNSIGNED NOT NULL DEFAULT '0'
     COMMENT 'size of source in bytes',

  INDEX(class_fk),
  FOREIGN KEY(class_fk) REFERENCES Parameter_class (Parameter_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(vote_fk) REFERENCES Votes (vote_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB
COMMENT='a row defines an instance of one of the parameter types in Parameter_class';


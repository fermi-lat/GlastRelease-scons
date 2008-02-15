create table Ancillary
(Ancillary_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
 class_fk INT UNSIGNED NOT NULL
   COMMENT 'refers to row in Ancillary_class',

  instrument  varchar(16) not null default 'LAT'
     COMMENT 'e.g. LAT, EM2  ',

  source                 VARCHAR(255) NOT NULL 
     COMMENT 'path to file in Moot archive, relative to root',

  source_fmt             VARCHAR(32) NOT NULL
     COMMENT 'could be, e.g. xml',

  GMT_create_time         DATETIME 
  COMMENT 'When file this row refers to was created',

  creation_time         DATETIME 
  COMMENT 'When this row was created',

  creator               VARCHAR(32) NOT NULL
  COMMENT 'Who or what created the row',

  quality               ENUM('PROD', 'DEV', 'TEST', 'SUPSED', 'INVALID') 
                       NOT NULL,

  description           VARCHAR(255) NOT NULL
     COMMENT 'a place to mention anything special about this ancillary 
            instance',

##  tower                 TINYINT DEFAULT '-1' NOT NULL
##     COMMENT 'Ancillary info may be by tower or by entire detector',

  status                 ENUM('STARTED', 'ABORT', 'CREATED', 'INVALID')
                            NOT NULL DEFAULT 'STARTED'
COMMENT 'refers to status of entry, in case file has to be copied',

  inspect_state VARCHAR(32) NOT NULL DEFAULT 'NOT_INSPECTED' 
      COMMENT 'typical values are NOT_INSPECTED, PASSED, FAILED..',

  INDEX(class_fk),
  FOREIGN KEY(class_fk) REFERENCES Ancillary_class (Ancillary_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB
COMMENT='a row defines an instance of one of the ancillary types in Ancillary_class';

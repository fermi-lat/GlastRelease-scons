create table FSW_inputs
 (FSW_input_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  class_fk               INT UNSIGNED NOT NULL 
     COMMENT 'refers to row in FSW_class',

  FSW_id                 INT UNSIGNED NOT NULL
     COMMENT 'id from FSW logical table when binary has been created', 

  instrument  varchar(16) not null default 'LAT'
      COMMENT 'e.g. LAT, EM2  ',

  fmx_db  varchar(32) not null default 'none'
      COMMENT 'recognized values are development, production, none  ',

  status                 VARCHAR(32) NOT NULL
     COMMENT 'may take on values including "new", "added", "invalid"',

  source                 VARCHAR(255) NOT NULL 
     COMMENT 'path to source file in Moot archive, relative to root',

  source_fmt             VARCHAR(32) NOT NULL
     COMMENT 'could be, e.g. xml',

  fmx_path                VARCHAR(255) NOT NULL DEFAULT 'unknown'
     COMMENT 'path within fmx archive relative to root; used in fof',

  creation_time         DATETIME
  COMMENT 'When this row was created',

  creator               VARCHAR(32) NOT NULL,

  description           VARCHAR(255) NOT NULL
     COMMENT 'keep track of precursors with list of keys  (Parameters or FSW_inputs) ',
  checksum              INT UNSIGNED  
     COMMENT 'provides quick way for mostly-reliable diff',
  size                  INT UNSIGNED NOT NULL DEFAULT '0'
     COMMENT 'size of source in bytes',
  INDEX(class_fk),
  FOREIGN KEY(class_fk) REFERENCES FSW_class (FSW_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB
COMMENT='a row defines an instance of one of the FSW_input types in FSW_class';


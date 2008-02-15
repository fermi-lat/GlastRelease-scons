create table FSW_class
 (FSW_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name              VARCHAR(32) NOT NULL UNIQUE,
  from_precursors       VARCHAR(32) NOT NULL
   COMMENT 'indicate how to create FSW_input from precursors',
  to_binary             VARCHAR(32) NOT NULL
   COMMENT 'indicate FSW service to invoke to create binary',
  description   VARCHAR(255),
  active_state VARCHAR(32) NOT NULL DEFAULT 'ACTIVE' 
    COMMENT 'typical values are ACTIVE, INACTIVE'
)  TYPE=InnoDB
COMMENT='Kinds of FSW inputs used by MOOT';

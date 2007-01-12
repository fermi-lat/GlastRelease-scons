create table FSW_to_Parameters
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  Parameter_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to row in Parameters',

  FSW_fk      INT UNSIGNED NOT NULL
    COMMENT 'refers to row in FSW_inputs',

  status      VARCHAR(32) not null default 'ACTIVE'
     COMMENT 'any value other than default means unusable',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
  INDEX(FSW_fk),
  FOREIGN KEY(FSW_fk) REFERENCES FSW_inputs (FSW_input_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(Parameter_fk) REFERENCES Parameters (parm_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate FSW_input to its Parameters';


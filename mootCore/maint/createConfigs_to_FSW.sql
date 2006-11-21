create table Configs_to_FSW
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  config_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to row in Configs',

  FSW_fk      INT UNSIGNED NOT NULL
    COMMENT 'refers to row in FSW_inputs',

  status      VARCHAR(32) not null default 'active'
     COMMENT 'any value other than default means unusable',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
  INDEX(config_fk),
  FOREIGN KEY(config_fk) REFERENCES Configs (config_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(FSW_fk) REFERENCES FSW_inputs (FSW_input_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Config to its FSW_inputs';


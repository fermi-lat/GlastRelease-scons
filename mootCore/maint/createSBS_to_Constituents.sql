create table SBS_to_Constituents
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  sbs_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to row in FSW_inputs',

  Constituents_fk      INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Constituents',

  status      VARCHAR(32) not null default 'active'
     COMMENT 'any value other than default means unusable',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
  INDEX(sbs_fk),
  FOREIGN KEY(sbs_fk) REFERENCES FSW_inputs (FSW_input_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(Constituents_fk) REFERENCES Constituents (prim_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate SBS key to its Constituents';


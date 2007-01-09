create table Configs_to_Parameters
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  config_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to row in Configs',
  Parameter_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to row in Parameters rqsted when config was built',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
  INDEX(config_fk),
  FOREIGN KEY(config_fk) REFERENCES Configs (config_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(Parameter_fk) REFERENCES Parameters (parm_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Config to requested (not nec. used) parameter files';


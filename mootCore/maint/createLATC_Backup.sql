create table LATC_Backup
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  FSW_fk   INT UNSIGNED NOT NULL 
    COMMENT 'refers to latc-type FSW class',

  pclass_fk      INT UNSIGNED
    COMMENT 'refers to row in Parameter_class for auto-generated backup parm',

  active_state       enum('ACTIVE', 'INACTIVE') NOT NULL default 'ACTIVE'
     COMMENT 'provide a way to disable the association',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
    FOREIGN KEY(FSW_fk) REFERENCES FSW_class(FSW_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(pclass_fk) REFERENCES Parameter_class (Parameter_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  UNIQUE INDEX (FSW_fk)
 ) TYPE=InnoDB COMMENT= 
'Relate latc component to corresponding backup source class';


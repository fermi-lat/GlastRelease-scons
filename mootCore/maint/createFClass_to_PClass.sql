create table FClass_to_PClass
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  PClass_fk INT UNSIGNED NOT NULL UNIQUE
    COMMENT 'refers to row in Parameter_class table',
  FClass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in FSW_class table',
  status    ENUM('ACTIVE','TEST','INACTIVE') NOT NULL DEFAULT 'TEST'
    COMMENT 'status of relation: should it be used to build new Config',
  one_to_one TINYINT  NOT NULL DEFAULT '1' 
    COMMENT 'value is 1 iff parameter file is complete source for FSW input',
  FOREIGN KEY(PClass_fk) REFERENCES Parameter_class (Parameter_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(FClass_fk) REFERENCES FSW_class (FSW_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Parameter_class and FSW_class';


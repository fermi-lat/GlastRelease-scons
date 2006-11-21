create table PClass_to_AClass
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
    COMMENT 'primary key',
  PClass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Parameter_class table',
  AClass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Ancillary table',
  status    ENUM('ACTIVE','TEST','INACTIVE') NOT NULL DEFAULT 'TEST'
    COMMENT 'status of relation: should it be used to build new Config',
  FOREIGN KEY(PClass_fk) REFERENCES Parameter_class (Parameter_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(AClass_fk) REFERENCES Ancillary_class (Ancillary_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Parameter_class and Ancillary_class';


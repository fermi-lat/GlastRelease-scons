create table OClass_to_AClass
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  AClass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Ancillary_class table',
  OClass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in OfflineCalib table',
  status    ENUM('ACTIVE','TEST','INACTIVE') NOT NULL DEFAULT 'TEST'
    COMMENT 'status of relation: should it be used to build new Config',
  INDEX(OClass_fk),
  FOREIGN KEY(OClass_fk) REFERENCES Offline_class (Offline_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(AClass_fk) REFERENCES Ancillary_class (Ancillary_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Offline_class and Ancillary_class';


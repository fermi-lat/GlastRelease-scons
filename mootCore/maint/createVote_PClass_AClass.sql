create table Vote_PClass_AClass
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
    COMMENT 'primary key',
  vote_fk   INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Votes table',
  pclass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Parameter_class table',
  aclass_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Ancillary_class table',
  a_alias   VARCHAR(32) NOT NULL
    COMMENT 'alias name for ancillary file',
  status    ENUM('ACTIVE','TEST','INACTIVE') NOT NULL DEFAULT 'TEST'
    COMMENT 'status of relation. Expect it will always be ok (ACTIVE)',
  creation_time DATETIME
    COMMENT 'time (GMT) of initial insert of this row',
  creator  VARCHAR(32) NOT NULL
    COMMENT 'OS username of person making the insert of this row',

  INDEX(vote_fk),
  FOREIGN KEY(vote_fk) REFERENCES Votes (vote_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(pclass_fk) REFERENCES Parameter_class (Parameter_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(aclass_fk) REFERENCES Ancillary_class (Ancillary_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate Parameter_class and Ancillary_class';


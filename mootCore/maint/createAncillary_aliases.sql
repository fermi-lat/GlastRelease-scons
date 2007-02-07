create table Ancillary_aliases
(key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
 COMMENT "primary key", 
 name  VARCHAR(32) NOT NULL
 COMMENT "things like 'nominal' ",
 aclass_fk INT UNSIGNED NOT NULL default '1'
  COMMENT 'refers to row in Ancillary_class table',
 tower TINYINT NOT NULL default '-1'
  COMMENT "use -1 for 'all' or 'not applicable',
 ancillary_fk INT UNSIGNED NOT NULL
  COMMENT 'points to ancillary file which is the resolution of the
  alias name for specified class, tower'
 last_modified_time DATETIME
 COMMENT 'last update. In practice, last time ancillary_fk was written',
 last_modifier       VARCHAR(32) NOT NULL
  COMMENT 'OS username of person who did last update',
 FOREIGN KEY(aclass_fk) REFERENCES Ancillary_class(Ancillary_class_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
 FOREIGN KEY(ancillary_fk) REFERENCES Ancillary(Ancillary_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
 UNIQUE INDEX (name, aclass_fk, tower)
) TYPE=InnoDB
COMMENT 'intent files always refer to ancillary files by alias;

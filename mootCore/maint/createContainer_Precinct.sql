create table Container_Precinct
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
    COMMENT 'primary key',
  ctn_fk   INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Votes table for container vote',
  precinct_fk INT UNSIGNED NOT NULL
    COMMENT 'refers to individual precinct mentioned in container',
  v_alias   VARCHAR(32) NOT NULL
    COMMENT 'alias name for individual precinct vote file',
  status    ENUM('ACTIVE','TEST','INACTIVE') NOT NULL DEFAULT 'TEST'
    COMMENT 'status of relation. Expect it will always be ok (ACTIVE)',
  creation_time DATETIME
    COMMENT 'time (GMT) of initial insert of this row',
  creator  VARCHAR(32) NOT NULL
    COMMENT 'OS username of person making the insert of this row',

  INDEX(ctn_fk),
  FOREIGN KEY(precinct_fk) REFERENCES Precincts (precinct_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
  FOREIGN KEY(ctn_fk) REFERENCES  Votes (vote_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate vote container file and individual votes';


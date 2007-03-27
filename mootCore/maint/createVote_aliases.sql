create table Vote_aliases
(vote_aliases_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
 COMMENT "primary key", 
 name  VARCHAR(32) NOT NULL
 COMMENT "things like 'nominal' ",
 name  VARCHAR(32) NOT NULL
 COMMENT "things like 'nominal' ",
 precinct_fk INT UNSIGNED NOT NULL default '1'
  COMMENT 'refers to row in Precincts table',
 vote_fk INT UNSIGNED NOT NULL
  COMMENT 'points to vote file which is the resolution of the
  alias name for specified class (i.e., precinct)',
 last_modified_time DATETIME
 COMMENT 'last update. In practice, last time ancillary_fk was written',
 last_modifier       VARCHAR(32) NOT NULL
  COMMENT 'OS username of person who did last update',
 FOREIGN KEY(precinct_fk) REFERENCES Precincts(precincts_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
 FOREIGN KEY(vote_fk) REFERENCES Votes(vote_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,
 UNIQUE INDEX (name, precinct_fk)
) TYPE=InnoDB
COMMENT 'global intent files may refer to precinct-wide intent files (votes) by alias';

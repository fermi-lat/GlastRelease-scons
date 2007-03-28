create table Votes
(vote_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
 COMMENT "primary key", 
 name  VARCHAR(32) NOT NULL
  COMMENT "need not be unique",
 precinct_fk INT UNSIGNED NOT NULL  default '1'
  COMMENT "precinct this vote pertains to",
 source                  VARCHAR(255) NOT NULL 
     COMMENT 'path to file in Moot archive, relative to root',
 status    VARCHAR(32) NOT NULL default 'OK'
  COMMENT 'set to some value other than OK to indicate it 
           should not be used',
description           VARCHAR(255) NOT NULL
  COMMENT 'a place to mention anything special about the Vote',
 creation_time DATETIME
  COMMENT 'time (GMT) of initial insert of this row',
 creator  VARCHAR(32) NOT NULL
  COMMENT 'OS username of person making the insert of this row',
 INDEX(precinct_fk),
 FOREIGN KEY(precinct_fk) REFERENCES Precincts (precinct_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
) TYPE=InnoDB


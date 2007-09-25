create table Votes
(vote_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
 COMMENT "primary key", 
 precinct_fk INT UNSIGNED NOT NULL  default '1'
  COMMENT "precinct this vote pertains to",
 source                  VARCHAR(255) NOT NULL 
     COMMENT 'path to file in Moot archive, relative to root',
  status                 ENUM('STARTED', 'ABORT', 'CREATED', 'INVALID')
                            NOT NULL DEFAULT 'STARTED'
  COMMENT 'refers to status of entry.  Not usable until it has been
           parsed, copied and has minimally sensible contents',
  instrument  varchar(16) not null default 'LAT'
     COMMENT 'e.g. LAT, EM2  ',

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


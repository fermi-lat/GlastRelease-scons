create table Precincts
(key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
 COMMENT "primary key", 
 name  VARCHAR(32) NOT NULL UNIQUE
 COMMENT "ACD_Mode, ACD_Bias, CAL_Timing, ...",
 delegate_fk INT UNSIGNED NOT NULL default '1'
  COMMENT 'refers to row in Delegate table',
 status    VARCHAR(32) NOT NULL default 'OK'
  COMMENT 'set to some value other than OK to indicate it 
           cannot be referred to by new entries in other tables',
 creation_time DATETIME
  COMMENT 'time (GMT) of initial insert of this row',
 creator  VARCHAR(32) NOT NULL
  COMMENT 'OS username of person making the insert of this row',
 INDEX(delegate_fk),
 FOREIGN KEY(delegate_fk) REFERENCES Delegates (delegate_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
) TYPE=InnoDB
COMMENT 'a row defines the smallest unit of configuration management';
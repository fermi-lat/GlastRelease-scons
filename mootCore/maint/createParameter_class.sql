##  Should get rid of delegate_fk one of these days.  Precinct_fk is the
## important thing

create table Parameter_class
(Parameter_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name              VARCHAR(32) NOT NULL UNIQUE,
  description       VARCHAR(255) NOT NULL,
  active_state      ENUM('ACTIVE','INACTIVE') NOT NULL default 'ACTIVE'
    COMMENT 'this field allows for evolution of collection of Parameter classes',
  precinct_fk INT UNSIGNED  COMMENT 
   'who manages this parameter class',
  FOREIGN KEY(precinct_fk) REFERENCES Precincts (precinct_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE

  delegate_fk INT UNSIGNED  COMMENT 
   'who manages this parameter class',
  FOREIGN KEY(delegate_fk) REFERENCES Delegates (delegate_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE

)  TYPE=InnoDB
COMMENT='kinds of parameter sets used by MOOT';


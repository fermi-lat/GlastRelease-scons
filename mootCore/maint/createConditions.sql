create table Conditions
(condition_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
  COMMENT "primary key",

 name       VARCHAR(255) NOT NULL 
  COMMENT "expected to be close to unique",



# Put a comment in the script

status    ENUM('OK', 'SUPSED', 'INVALID') NOT NULL
  COMMENT 'is this conditions description usable or not',

instrument  varchar(16) NOT NULL default 'LAT'
   COMMENT 'e.g. LAT, EM2  ',

description           VARCHAR(255) NOT NULL
  COMMENT 'a place to mention anything special about this CondSet',

delegate_fk INT UNSIGNED NOT NULL COMMENT 
'a non-trivial Condition set applies to one delegate; convention
 for delegates is that 1 stands for -none- ',

# Delegate-specific fields

legain    ENUM('0', '1', '2', '3', '4', '5', '6', '7', 'none') NOT NULL
    default 'none' 
 COMMENT 'CAL low-energy gain setting, normally 5',

hegain   
ENUM('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15',
     'none') NOT NULL default 'none' 
   COMMENT 'CAL high-energy gain setting, normally 15 for flight',

temperature FLOAT  NOT NULL default '0' COMMENT 'degrees C',

temperature_tol FLOAT NOT NULL default '0' 
   COMMENT 'degrees C; if 0 temperature is unused' ,

creation_time   DATETIME 
  COMMENT 'When this row was created',

creator   VARCHAR(32) NOT NULL
  COMMENT 'username of creator',

  FOREIGN KEY(delegate_fk) REFERENCES Delegates (delegate_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE

  ) TYPE=InnoDB
COMMENT='a row describes conditions for parameters belonging to a particular delegate-mode';


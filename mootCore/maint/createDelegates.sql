create table Delegates
(delegate_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
  COMMENT "primary key",

 name       VARCHAR(16) NOT NULL UNIQUE
  COMMENT "CAL, TKR, etc.",

# Put a comment in the script

status    ENUM('OK',  'INACTIVE') NOT NULL
  COMMENT 'should be OK for a real delegate',

description           VARCHAR(255) NOT NULL
  COMMENT 'a place to mention anything special about the delegate',


creation_time   DATETIME 
  COMMENT 'When this row was created',

creator   VARCHAR(32) NOT NULL
  COMMENT 'username of creator'


  ) TYPE=InnoDB
COMMENT='a row describes a management entity, like CAL';



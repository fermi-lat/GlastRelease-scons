create table Ancillary_class
(Ancillary_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                     VARCHAR(32) NOT NULL UNIQUE,
  description       VARCHAR(255) NOT NULL
)  TYPE=InnoDB
COMMENT='ancillary file classes';

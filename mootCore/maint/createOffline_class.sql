create table Offline_class
(Offline_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name              VARCHAR(32) NOT NULL UNIQUE,
  config            TINYINT     NOT NULL DEFAULT '0'
    COMMENT 'True (=1) if calib is set-up data; false if it is response info',
  description       VARCHAR(255) NOT NULL
) TYPE=InnoDB
COMMENT='Kinds of offline calibrations';

# Contains acq. information which may vary depending on mode.
# At a minimum this includes filter configuration
# Would like to have another column indicating..   
# [?? too bad I didn't finish the above thought]

# Might change name FSW_input_fk  --> FSW_fk (for consistency)
#  or might change it to refer to new Constituents table rather than FSW_inputs
#  table
create table Acq_mode
 (acq_mode_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY
    COMMENT "primary key",
  acq_summary_fk INT UNSIGNED NOT NULL
    COMMENT "identifies acquisition",
  fmx_logical INT UNSIGNED NOT NULL
    COMMENT "logical key for config file",
  FSW_input_fk      INT UNSIGNED 
    COMMENT "refers to row in FSW_inputs with FSW_inputs.FSW_id = our fmx_logical, if found",
  mode       VARCHAR(32) NOT NULL
    COMMENT "e.g TOO, ARR, normal acq",
  INDEX(acq_summary_fk),
  UNIQUE INDEX(acq_summary_fk, fmx_logical, mode),
  FOREIGN KEY(FSW_input_fk) REFERENCES FSW_inputs(FSW_input_key)
    ON DELETE SET NULL
    ON UPDATE CASCADE,
  FOREIGN KEY(acq_summary_fk) REFERENCES Acq_summary(acq_summary_key)
    ON DELETE CASCADE
    ON UPDATE CASCADE
 ) TYPE=InnoDB
COMMENT='a row for each unique triple (acq_summary_fk, fmx_logical, mode)';

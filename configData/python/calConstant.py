

"""
Global calorimeter constant values
"""


__facility__      = "Online"
__abstract__      = "Global calorimeter constant values"
__author__        = "D.L. Wood"
__date__          = "$Date$"
__version__       = "$Revision$, $Author$"
__release__       = "$Name$"
__credits__       = "NRL code 7650"


# Number of possible calorimeters TEM's
NUM_TEM = 16

# Number of cal cable controllers
NUM_GCCC = 4


# Number of cal row controllers
NUM_GCRC = 4

# Number of calorimeter FE's
NUM_GCFE = 12

# Total number of GCCC in cal
GLOBAL_N_GCCC = NUM_TEM*NUM_GCCC

# Total Number of GCRC in cal
GLOBAL_N_GCRC = GLOBAL_N_GCCC*NUM_GCRC

# Total Number of GCFE's in Cal
GLOBAL_N_GCFE = GLOBAL_N_GCRC*NUM_GCFE

# nominal threhsold dac slopes
NOMINAL_THRESH_DAC_SLOPE = {}
NOMINAL_THRESH_DAC_SLOPE["log_acpt"]     = .33
NOMINAL_THRESH_DAC_SLOPE["fle_dac"]      = 4
NOMINAL_THRESH_DAC_SLOPE["fhe_dac"]      = 50
NOMINAL_THRESH_DAC_SLOPE["rng_uld_dac"]  = 28

# units for THRESH_DAC_SLOPE
THRESH_DAC_SLOPE_UNITS = {}
THRESH_DAC_SLOPE_UNITS["log_acpt"]    = "MeV/DAC"
THRESH_DAC_SLOPE_UNITS["fle_dac"]     = "MeV/DAC"
THRESH_DAC_SLOPE_UNITS["fhe_dac"]     = "MeV/DAC"
THRESH_DAC_SLOPE_UNITS["rng_uld_dac"] = "ADC/DAC"


# main register name associated with cal threshold precincts
THRESH_REG_NAME = {}
THRESH_REG_NAME["CAL_LAC"] = "log_acpt"
THRESH_REG_NAME["CAL_FLE"] = "fle_dac"
THRESH_REG_NAME["CAL_FHE"] = "fhe_dac"
THRESH_REG_NAME["CAL_ULD"] = "rng_uld_dac"


def getCCCId(index):
    """
    return tuple (tem,ccc) of cal component indices for given CCC flat index
    """
    ccc = index % NUM_GCCC
    tem = int(index/NUM_GCCC)

    return (tem,ccc)

def getCRCId(index):
    """
    return tuple (tem,ccc,crc) of cal component indices for given CRC flat index
    """
    crc = index % NUM_GCRC
    index /= NUM_GCRC
    (tem,ccc) = getCCCId(index)

    return (tem,ccc,crc)

def getCFEId(index):
    """
    return tuple (tem,ccc,crc,cfe) of cal component indices for given CFE flat index
    """
    cfe = index % NUM_GCFE
    index /= NUM_GCFE
    (tem,ccc,crc) = getCRCId(index)

    return (tem,ccc,crc,cfe)



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

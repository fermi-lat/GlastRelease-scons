"""
Classes to represent CAL DAC settings XML documents.
"""


__facility__  = "Offline"
__abstract__  = "Classes to represent CAL DAC settings XML documents"
__author__    = "D.L.Wood"
__date__      = "$Date$"
__version__   = "$Revision$, $Author$"
__release__   = "$Name$"
__credits__   = "NRL code 7650"



import logging
import xml.dom.minidom

import Numeric

import calXML
from calExcept import *



MODE_CREATE     = calXML.MODE_CREATE
MODE_READONLY   = calXML.MODE_READONLY


DAC_FILE_TYPES = ('fle_dac', 'fhe_dac', 'log_acpt', 'rng_uld_dac')
ENG_FILE_TYPES = ('thrBias', 'adc2nrg')




class calSnapshotXML(calXML.calXML):
    """
    CAL hardware configuration snapshot XML data file class.
    
    This class provides methods for accessing CAL HW configuration
    data stored in XML snapshot file format.
    """

    def __init__(self, fileName, mode = MODE_READONLY):
    
        calXML.calXML.__init__(self, fileName, mode)

        
    def read(self, dacName):
        """
        Read CAL data from a snapshot XML file

        Param: dacName - The name of the DAC element to read.  The following
        HW items are supported:

            fle_dac - CAL FLE DAC setttings
            fhe_dac - CAL FHE DAC settings
            log_acpt - CAL LAC DAC settings
            rng_uld_dac - CAL ULD DAC settings        

        Returns: A Numeric array of data (16, 8, 2, 12) read from the
        XML file.
        """

        # create empty DAC data array

        dacData = Numeric.zeros((16, 8, 2, 12), Numeric.Int16)

        # find <fle_dac> elements

        doc = self.getDoc()

        latList = doc.getElementsByTagName('GLAT')
        latLen = len(latList)
        if latLen != 1:
            raise calFileReadExcept, "found %d <GLAT> elements, expected 1" % latLen
        temList = latList[0].getElementsByTagName('GTEM')
        temLen = len(temList)
        if temLen > 16:
            raise calFileReadExcept, "found %d <GTEM> elements, expected <= 16" % temLen
        for t in temList:
            tem = int(t.getAttribute('ID'))
            if tem < 0 or tem > 16:
                raise calFileReadExcept, "<GTEM> ID attribute value %d, expected (0 - 15)" % tem
            cccList = t.getElementsByTagName('GCCC')
            cccLen = len(cccList)
            if cccLen > 4:
                raise calFileReadExcept, "found %d <GCCC> elements, expected <= 4" % cccLen
            for c in cccList:
                ccc = int(c.getAttribute('ID'))
                if ccc < 0 or ccc > 3:
                    raise calFileReadExcept, "<GCCC> ID attribute value %d, expected (0 - 3)" % ccc
                rcList = c.getElementsByTagName('GCRC')
                rcLen = len(rcList)
                if rcLen > 4:
                    raise calFileReadExcept, "found %d <GCRC> elements, expected <= 4" % rcLen
                for r in rcList:
                    rc = int(r.getAttribute('ID'))
                    if rc < 0 or rc > 3:
                        raise calFileReadExcept, "<GCRC> ID attribute value %d, expected (0 - 3)" % rc
                    feList = r.getElementsByTagName('GCFE')
                    feLen = len(feList)
                    if feLen > 12:
                        raise calFileReadExcept, "found %d <GCFE> elements, expected <= 12" % feLen
                    for f in feList:
                        fe = int(f.getAttribute('ID'))
                        if fe < 0 or fe > 11:
                            raise calFileReadExcept, "<GCFE> ID attribute value %d, expected (0 - 11)" % fe
                        dacList = f.getElementsByTagName(dacName)
                        dacLen = len(dacList)
                        if dacLen != 1:
                            raise calFileReadExcept, "found %d %s elements, expected 1" % (dacLen, dacName)
                        d = dacList[0]
                        dd = d.childNodes[0]
                        dac = int(dd.data.strip(), 16)
                        if (ccc % 2) != 0:
                            row = (rc + 4)
                        else:
                            row = rc
                        if ccc < 2:
                            end = 1
                        else:
                            end = 0
                        dacData[tem, row, end, fe] = dac

        return dacData

    

class calDacXML(calSnapshotXML):   
    """
    CAL DAC configuration XML data file class.
    
    This class provides methods for accessing CAL DAC configuration
    data stored in XML format.  These files are fragments of
    full HW snapshot records.

    The following file types are supported:

        fle_dac - CAL FLE DAC setttings
        fhe_dac - CAL FHE DAC settings
        log_acpt - CAL LAC DAC settings
        rng_uld_dac - CAL ULD DAC settings
        thrBias - CAL threshold bias values
        adc2nrg - CAL ADC to energy conversion values
    """

    def __init__(self, fileName, dacName, mode = MODE_READONLY):
        """
        Open a CAL DAC configuration XML file

        \param fileName The XML file name.
        \param dacName The name of the bottom level XML data element.
        \param mode The file access mode (MODE_READONLY or MODE_CREATE).
        """

        calSnapshotXML.__init__(self, fileName, mode)  

        if dacName not in DAC_FILE_TYPES:
            raise calFileOpenExcept, "DAC name %s not supported" % str(dacName)
        self.__dacName = dacName


    def read(self):
        """
        Read data from a CAL DAC configuration XML file

        Returns: A Numeric array of data (16, 8, 2, 12) read from the
        XML file.
        """

        return calSnapshotXML.read(self, self.__dacName)


    def write(self, dacData, lrefgain = None, hrefgain = None, tems = (0,)):
        """
        Write CAL data to a snapshot fragment XML file

        Param: dacData - A Numeric array of data (16, 8, 2, 12) to write to the
        XML file.
        Param: hrefgain 
        Param: tems - A list of TEM ID values to include in the output data set.
        """

        doc = self.getDoc()

        # insert <LATdoc> element

        ld = doc.createElement('LATdoc')
        ld.setAttribute('name', self.getName())
        doc.appendChild(ld)

        # insert <configuration> element

        ce = doc.createElement('configuration')
        ce.setAttribute('name', '')
        s = '[\'GCCC\',\'GCRC\',\'GCFE\',\'%s\']' % self.__dacName  
        ce.setAttribute('hierarchy', s)
        ce.setAttribute('type', 's')
        ce.setAttribute('shape', '(16,8,2,12)')
        ce.setAttribute('version', 'NA')
        if lrefgain is not None:
            ce.setAttribute('leRefGain', str(lrefgain))
        if hrefgain is not None:
            ce.setAttribute('heRefGain', str(hrefgain))
        ld.appendChild(ce)

        # insert <GLAT> element  
            
        gl = doc.createElement('GLAT')
        ld.appendChild(gl)

        for tem in tems:
                        
            # insert <GTEM> elements

            gt = doc.createElement('GTEM')
            gt.setAttribute('ID', str(tem))
            gl.appendChild(gt)
            
            for ccc in range(4):

                # insert <GCCC> elements

                gc = doc.createElement('GCCC')
                gc.setAttribute('ID', str(ccc))
                gt.appendChild(gc)                

                for rc in range(4):

                    # insert <GCRC> elements

                    gr = doc.createElement('GCRC')
                    gr.setAttribute('ID', str(rc))
                    gc.appendChild(gr)
                
                    # translate index

                    if (ccc % 2) != 0:
                        row = (rc + 4)
                    else:
                        row = rc
                    if ccc < 2:
                        end = 1
                    else:
                        end = 0                                            

                    for fe in range(12):

                        # insert <GCFE> elements

                        gf = doc.createElement('GCFE')
                        gf.setAttribute('ID', str(fe))
                        gr.appendChild(gf)

                        # insert <xxx> elements

                        dv = doc.createElement(self.__dacName)
                        t = doc.createTextNode('0x%x' % int(dacData[tem, row, end, fe]))
                        dv.appendChild(t)
                        gf.appendChild(dv)
          
        
        # write output XML file

        self.writeFile()



class calEnergyXML(calXML.calXML):
    """
    CAL ADC to energy conversion data.
    
    This class provides methods for accessing CAL ADC to energy
    data stored in XML pseudo-snapshot file format.
    """

    def __init__(self, fileName, engName, mode = MODE_READONLY):
        """
        Open a CAL DAC configuration XML file

        \param fileName The XML file name.
        \param engName The name of the bottom level XML data element.
        \param mode The file access mode (MODE_READONLY or MODE_CREATE).
        """

        calXML.calXML.__init__(self, fileName, mode)  

        if engName not in ENG_FILE_TYPES:
            raise calFileOpenExcept, "ENG name %s not supported" % str(engName)
        self.__engName = engName

        
    def read(self):
        """
        Read CAL data from a snapshot XML file

        Returns: A Numeric array of data (16, 8, 2, 12, 2) read from the
        XML file.  The last dimension is as follows:
            0 = LE conversion value
            1 = HE conversion value
        """

        doc = self.getDoc()        

        # create empty energy data array

        engData = Numeric.zeros((16, 8, 2, 12, 2), Numeric.Float32)

        # find elements

        latList = doc.getElementsByTagName('GLAT')
        latLen = len(latList)
        if latLen != 1:
            raise calFileReadExcept, "found %d <GLAT> elements, expected 1" % latLen
        temList = latList[0].getElementsByTagName('GTEM')
        temLen = len(temList)
        if temLen > 16:
            raise calFileReadExcept, "found %d <GTEM> elements, expected <= 16" % temLen
        for t in temList:
            tem = int(t.getAttribute('ID'))
            if tem < 0 or tem > 16:
                raise calFileReadExcept, "<GTEM> ID attribute value %d, expected (0 - 15)" % tem
            if self.__engName == 'adc2nrg':
                eName = 'low_hi_nrg'
            else:
                eName = 'trhBias'
            eList = t.getElementsByTagName(eName)
            eLen = len(eList)
            if eLen > 2:
                raise calFileReadExcept, "found %d <GCCC> elements, expected <= 2" % eLen
            for e in eList:
                eng = int(e.getAttribute('ID'))
                if eng < 0 or eng > 1:
                    raise calFileReadExcept, "%s ID attribute value %d, expected (0 - 1)" % (eName, eng)
                cccList = e.getElementsByTagName('GCCC')
                cccLen = len(cccList)
                if cccLen > 4:
                    raise calFileReadExcept, "found %d <GCCC> elements, expected <= 4" % cccLen
                for c in cccList:
                    ccc = int(c.getAttribute('ID'))
                    if ccc < 0 or ccc > 3:
                        raise calFileReadExcept, "<GCCC> ID attribute value %d, expected (0 - 3)" % ccc
                    rcList = c.getElementsByTagName('GCRC')
                    rcLen = len(rcList)
                    if rcLen > 4:
                        raise calFileReadExcept, "found %d <GCRC> elements, expected <= 4" % rcLen
                    for r in rcList:
                        rc = int(r.getAttribute('ID'))
                        if rc < 0 or rc > 3:
                            raise calFileReadExcept, "<GCRC> ID attribute value %d, expected (0 - 3)" % rc
                        feList = r.getElementsByTagName('GCFE')
                        feLen = len(feList)
                        if feLen > 12:
                            raise calFileReadExcept, "found %d <GCFE> elements, expected <= 12" % feLen
                        for f in feList:
                            fe = int(f.getAttribute('ID'))
                            if fe < 0 or fe > 11:
                                raise calFileReadExcept, "<GCFE> ID attribute value %d, expected (0 - 11)" % fe
                            dacList = f.getElementsByTagName(self.__engName)
                            dacLen = len(dacList)
                            if dacLen != 1:
                                raise calFileReadExcept, "found %d %s elements, expected 1" % (dacLen, self.__engName)
                            d = dacList[0]
                            dd = d.childNodes[0]
                            dac = float(dd.data.strip())
                            if (ccc % 2) != 0:
                                row = (rc + 4)
                            else:
                                row = rc
                            if ccc < 2:
                                end = 1
                            else:
                                end = 0
                            engData[tem, row, end, fe, eng] = dac

        return engData


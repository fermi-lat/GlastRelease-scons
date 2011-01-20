import pSafeLogger
logger = pSafeLogger.getLogger('pXmlWriter')

import sys
import os
import time

from pAsciiWriter import pAsciiWriter

class pXmlWriter(pAsciiWriter):

    
    PROLOG = '<?xml version="1.0" ?>'

    def __init__(self, outputFilePath = None):
        pAsciiWriter.__init__(self, outputFilePath)

    def openFile(self, outputFilePath, mode = 'w'):
        pAsciiWriter.openFile(self, outputFilePath, mode)
        self.writeLine(self.PROLOG)
        self.newLine()

    def writeComment(self, line):
        self.writeLine('<!-- %s -->' % line)

    def getOpenTagBlock(self, tagName, attributesDict = {}, close = False):
        block = '<%s' % tagName
        for (key, value) in attributesDict.items():
            block += ' %s="%s"' % (key, value)
        if not close:
            block += '>'
        else:
            block += '/>'
        return block

    def getCloseTagBlock(self, tagName):
        return '</%s>' % tagName

 
    def getTagBlock(self, tagName, attributesDict = {}, value = None):
        if value is None:
            return self.getOpenTagBlock(tagName, attributesDict, True)
        else:
            return '%s%s%s' % (self.getOpenTagBlock(tagName, attributesDict),\
                               value, self.getCloseTagBlock(tagName))

    def openTag(self, tagName, attributesDict = {}, close = False):
        self.writeLine(self.getOpenTagBlock(tagName, attributesDict, close))

        
    def closeTag(self, tagName):
        self.writeLine(self.getCloseTagBlock(tagName))


    def writeTag(self, tagName, attributesDict, value = None):
        self.writeLine(self.getTagBlock(tagName, attributesDict, value))


    def writeBlock(self,attributesDict = {}):
        block = '<' 
        for (key, value) in attributesDict.items():
            block += '%s="%s"' % (key, value)
       
        block += '/>'
        return block

    def writeInfo(self,attributesDict):
        self.writeLine(self.writeBlock(attributesDict))

if __name__ == '__main__':
    writer = pXmlWriter('xml_test.xml')
    writer.openTag('alarmSummary')
    writer.indent()
    writer.writeTag('parameter', {'name': 'test', 'value': 3})
    writer.writeTag('output', {}, 10)
    writer.backup()
    writer.closeTag('alarmSummary')
    writer.closeFile()
